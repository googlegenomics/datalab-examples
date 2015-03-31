
# Navigating read and variant data

In this notebook we'll cover how retrieve read and variant data from the [Google
Genomics APIs](https://cloud.google.com/genomics/v1beta2/reference/).

This notebook is based on the ["Getting started in
python"](https://github.com/googlegenomics/getting-started/tree/master/python)
example.

### Setup

Before we get started, let's make sure we have an API key for making
authenticated requests.

*TODO(bryantd): figure out a simple auth story for these notebooks and use it
consistently throughout each*

*TODO(bryantd): cross-link/reference getting started with API notebook here*


    # TODO(bryantd): Come up with a better scheme for users to add their own
    # API key that avoids saving it in the notebook directly. It would be nice
    # to have the VM request an API key on startup actually and auto-set this
    # environment variable so that users never actually see it.
    
    API_KEY = os.environ.get('DNA_API_KEY', None)
    
    # Make sure that an API key was defined
    # todo: put instructions here for getting an api key
    if API_KEY is None:
        print 'Missing API key environment variable'

Now let's construct a Genomics API client object for making requests.


    from apiclient.discovery import build
    genomics = build('genomics', 'v1beta2', developerKey=API_KEY)

This `genomics` Python object will make HTTP requests on our behalf to the
Genomics API whenever we call the various resource methods.

Now we're ready to start fetching Genomics resources.

### Dataset selection

In this notebook, we'll examine a specific sample (NA12872) within the 1000
Genomes dataset. Let's define the set of identifiers we'll need to look up these
resources.


    dataset_id = '10473108253681171589' # This is the 1000 Genomes dataset ID
    sample = 'NA12872'
    reference_name = '22'
    reference_position = 51003835

### Find read group set ID

To begin, we need to lookup the ID of the read group set for our sample of
interest. We can use the `search` API of the `readgroupsets` resource to do
this. [Full documentation for the read group sets API can be found
here](https://cloud.google.com/genomics/v1beta2/reference/readgroupsets).

First, we'll create an appropriate search request:


    request = genomics.readgroupsets().search(
      # Our search criteria
      body={
        'datasetIds': [dataset_id],
        'name': sample
      },
      # The Set of fields we'd like returned
      fields='readGroupSets(id)')

And now we'll call the `execute()` method of our request, which will send an
HTTP request behind the scenes to execute our query against the Genomics API.


    response = request.execute()
    print response

    {u'readGroupSets': [{u'id': u'CMvnhpKTFhDu5oWzv-LJtIEB'}]}


We can now retrieve the ID that we're after from the response object:


    read_group_set_id = response['readGroupSets'][0]['id']
    print 'Read group set ID for sample %s is "%s"' % (sample, read_group_set_id)

    Read group set ID for sample NA12872 is "CMvnhpKTFhDu5oWzv-LJtIEB"


### Fetching reads

Now that we have the read group set ID, we can lookup the reads at the position
we are interested in.

We can do this by making another request to the Genomics API, following the same
steps as retrieving the read group set ID, except we'll access the **Reads**
resource this time.


    request = genomics.reads().search(
      body={'readGroupSetIds': [read_group_set_id],
            'referenceName': reference_name,
            'start': reference_position,
            'end': reference_position + 1,
            'pageSize': 1024},
      fields='alignments(alignment,alignedSequence)')
    
    reads = request.execute().get('alignments', [])
    
    # Note: This is simplistic - the cigar should be considered for real code
    bases = [read['alignedSequence'][
               reference_position - int(read['alignment']['position']['position'])]
             for read in reads]
    
    from collections import Counter
    print '%s bases on chromosome %s at %d are' % (sample, reference_name, reference_position)
    for base, count in Counter(bases).items():
      print '%s: %s' % (base, count)


    NA12872 bases on chromosome 22 at 51003835 are
    C: 1
    G: 13


### Select a Call Set ID

Now let's find the call set ID for our sample. We'll make a request to the
`Callsets` API this time.


    request = genomics.callsets().search(
      body={'variantSetIds': [dataset_id], 'name': sample},
      fields='callSets(id)')
    
    response = request.execute()
    
    call_set_id = response['callSets'][0]['id']


    print 'Found call set ID "%s"' % call_set_id

    Found call set ID "10473108253681171589-26"


### Fetching variants

Now that we have a call set ID, we can fetch the variants that overlap the
position we are interested in.


    request = genomics.variants().search(
      body={'callSetIds': [call_set_id],
            'referenceName': reference_name,
            'start': reference_position,
            'end': reference_position + 1},
      fields='variants(names,referenceBases,alternateBases,calls(genotype))')
    
    response = request.execute()
    
    # We'll look closer at the first variant
    variant = response['variants'][0]
    
    # Print the variant data we fetched
    for item in variant.items():
        print '%s => %s' % (item)

    calls => [{u'genotype': [1, 1]}]
    names => [u'rs131767', u'rs131767']
    alternateBases => [u'G']
    referenceBases => A


Finally, we can see what the called genotype of the variant is.

To map the genotype indices to bases, we need to concatenate the lists of
reference and alternate bases. For this case, this will result in `bases = ['A',
'G']`, so a genotype at `index=1` indicates `bases[1]` which would be G.


    variant_name = variants[0]['names'][0]
    bases = [variant['referenceBases']] + variant['alternateBases']
    print 'Bases: ', bases
    
    genotype = [bases[genotype_index] for genotype_index in variant['calls'][0]['genotype']]

    Bases:  [u'A', u'G']


And now we have the genotype for our variant


    print 'The called genotype is %s for %s' % (','.join(genotype), variant_name)

    The called genotype is G,G for rs131767



    

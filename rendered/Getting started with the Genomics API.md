
# Getting started with the Genomics APIs

In this notebook we'll cover how to make authenticated requests to the [Google
Genomics APIs](https://cloud.google.com/genomics/v1beta2/reference/), which will
include:
* Enabling the Genomics API for your project
* Obtaining an API key for making authenticated requests
* Creating a Python API client for issuing requests

We'll also cover how to make direct HTTP requests to the Genomics API.

## Setup

Before we can get started, we need to ensure that our VM is fully configured to
work against the Genomics APIs.

*TODO(bryantd): Move as much of this setup as possible into the Docker container
build.*

### Install Python libraries

We'll be using the [Google Python API client](https://github.com/google/google-
api-python-client) for interacting with Genomics APIs. We can install this
library, or any other 3rd-party Python libraries from the [Python Package Index
(PyPI)](https://pypi.python.org/pypi) using the `pip` package manager.

There are [50+ Google APIs](http://api-python-client-doc.appspot.com/) that you
can work against with the Google Python API Client, but we'll focus on the
Genomics API in this notebook.


    !pip install --upgrade google-api-python-client


### Get an API key

In order for us to work against public or private Genomics datasets, we'll need
to [get an API key](https://developers.google.com/api-client-
library/python/start/get_started#auth), which will authorize our requests.


    # TODO(bryantd): Come up with a better scheme for users to add their own
    # API key that avoids saving it in the notebook directly. It would be nice
    # to have the VM request an API key on startup actually and auto-set this
    # environment variable so that users never actually see it.

    API_KEY = os.environ.get('DNA_API_KEY', None)

    # Make sure that an API key was defined
    # todo: put instructions here for getting an api key
    if API_KEY is None:
        print 'Missing API key environment variable'

### Enable the Genomics API

You will need enable the Genomics API for your project if you have not done so
previously. The API section of the [Cloud Developer
Console](developers.google.com/console) will allow you to find and then enable
the Genomics API.

Once you've enabled the API, the next step is to construct a Python object that
we can use to make requests. The following snippet shows how we can create a
client for the Genomics API.


    from apiclient.discovery import build
    genomics = build('genomics', 'v1beta2', developerKey=API_KEY)

Now that we have a Python client for the Genomics API, we can access a variety
of different resources. For details about each available resource, see the [API
docs here](https://google-api-client-
libraries.appspot.com/documentation/genomics/v1beta2/python/latest/index.html).

Using our `genomics` client, we'll demonstrate fetching a Dataset resource by ID
(the 1000 genomes dataset in this case).

First, we need to construct a request object.


    request = genomics.datasets().get(datasetId='10473108253681171589')

Next, we'll send this request to the Genomics API by calling the
`request.execute()` method.


    response = request.execute()

The response object returned is simply a Python dictionary. Let's take a look at
the properties returned in the response.


    for entry in response.items():
        print "%s => %s" % entry

    isPublic => True
    id => 10473108253681171589
    name => 1000 Genomes
    projectNumber => 761052378059


Success! We can see the name of the specified Dataset and a few other pieces of
metadata.

Accessing other Genomics API resources will follow this same set of steps. The
full [list of available resources within the API is here](https://google-api-
client-
libraries.appspot.com/documentation/genomics/v1beta2/python/latest/index.html).
Each resource has details about the different verbs that can be applied (e.g.,
[Dataset methods](https://google-api-client-libraries.appspot.com/documentation/
genomics/v1beta2/python/latest/genomics_v1beta2.datasets.html)).

## Making direct HTTP requests to the API

It's also possible to directly access the Genomics API by making appropriately
formed HTTP requests. We'll use the `requests` library to make this easy.

Note that an API key is still used in this case, and is passed as a query string
parameter when making the HTTP request.


    import requests
    response = requests.get('https://www.googleapis.com/genomics/v1beta2/datasets/10473108253681171589',
                     params=dict(key=API_KEY))

The response object contains our requested data in the form of a JSON-encoded
string stored within the body of the response. We can transform this response to
a Python dictionary like so:


    import json
    data = json.loads(response.content)
    data




    {u'id': u'10473108253681171589',
     u'isPublic': True,
     u'name': u'1000 Genomes',
     u'projectNumber': u'761052378059'}



And access fields within the response data as expected:


    print data['name']

    1000 Genomes





#!/usr/bin/env python

import fnmatch
import os

from IPython.nbformat import current as nbformat
from IPython.nbconvert.exporters.markdown import MarkdownExporter

exporter = MarkdownExporter()

def get_notebook_filepaths(root_filepath):
  # Find all ipython notebooks within the notebooks folder
  notebooks = []
  for root, dirnames, filenames in os.walk(root_filepath):
    for filename in fnmatch.filter(filenames, '*.ipynb'):
        notebooks.append(os.path.join(root, filename))
  return notebooks


def convert_notebook_to_markdown(notebook_filepath, output_dirpath):
  with open(notebook_filepath) as f:
      nb = nbformat.reads_json(f.read())
  source, meta = exporter.from_notebook_node(nb)

  # TODO(bryantd): the output files for a notebook are not being generated at the moment
  # Issue is that the "ExtractOutputPreprocessor" needs to be run over the notebook and isn't
  # at the moment in this function. Maybe it's simpler to just shell out to nbconvert here
  # with the correct set of args at the moment, or use the NbConvertApp top-level ipy app object
  # to do the same without invoking the shell
  # See also: http://stackoverflow.com/questions/22037910

  # The GitHub Markdown rendering chokes if a newline doesn't exist between divs and tables
  # which causes the remainder of the Markdown document after such an instance to improperly render
  formatted_source = source.replace('</div><table>', '</div>\n<table>')

  nb_dirpath, nb_filepath = os.path.split(notebook_filepath)
  nb_filename, extension = os.path.splitext(nb_filepath)
  markdown_filepath = os.path.join(output_dirpath, nb_filename + '.md')

  with open(markdown_filepath, 'w+') as f:
      f.writelines(formatted_source)


if __name__ == '__main__':
  print "Converting IPython notebooks (.ipynb) to rendered Markdown..."
  notebook_filepaths = get_notebook_filepaths('notebooks')
  print 'Found %d notebooks to render' % len(notebook_filepaths)

  for nbpath in notebook_filepaths:
    print 'Converting .ipynb -> .md: ', nbpath
    convert_notebook_to_markdown(nbpath, 'rendered')

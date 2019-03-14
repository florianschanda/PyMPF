#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# PyMPF documentation build configuration file

import os
import sys
sys.path.insert(1, os.path.abspath('..'))

import mpf_version

# -- General configuration ------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.githubpages',
]

templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

project = 'PyMPF'
copyright = '2019, Florian Schanda'
author = 'Florian Schanda'

version = mpf_version.version
release = mpf_version.version

language = None

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#
# add_module_names = True

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True


# -- Options for HTML output ----------------------------------------------

html_theme = 'classic'
# html_theme_path = ['../sphinx']
# html_theme_options = {'collapsiblesidebar': True}

html_static_path = ['_static']

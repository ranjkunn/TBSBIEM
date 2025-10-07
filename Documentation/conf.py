# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'TBSBIEM'
copyright = '2025, Yenike Sharath Chandra Mouli and Ranjith Kunnath'
author = 'Yenike Sharath Chandra Mouli and Ranjith Kunnath'
release = '2025.0.0.0'

numfig = True
extensions = [
    'sphinx.ext.mathjax',  # For MathJax support
    # other extensions...
]


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

"""Public documentation only; import the installed package, never research scripts."""
import os
from importlib.metadata import version as package_version

project = "BReMAP"
release = package_version("bremap")
version = release
language = "en"
extensions = ["sphinx.ext.autodoc", "sphinx.ext.napoleon"]
autodoc_typehints = "none"
html_theme = "sphinx_rtd_theme"
html_title = f"BReMAP {release} documentation"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")
html_show_copyright = False
exclude_patterns = ["README.md"]
master_doc = "index"

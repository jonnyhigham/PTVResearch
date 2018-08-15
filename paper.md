---
title: 'PTVResearch: Robust particle tracking for all types of research applications'
tags:
  - MatLab
  - AppDesigner
  - Particle Tracking Velocimetry
  - Particle Image Velocimetry
  - Optical Flow
  - Lucas Kanade 
  - PODDEM
  - Proper Orthogonal Decomposition (POD)
authors:
  - name: Jonathan E. Higham
    orcid: 0000-0001-7577-0913
affiliations:
 - name: National Engergy Technology Labs., Department of Energy, Morgantown, WV.
date: 16 August 2018
bibliography: paper.bib
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

``Gala`` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for ``Gala`` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. ``Gala`` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the ``Astropy`` package [@astropy] (``astropy.units`` and
``astropy.coordinates``).

``Gala`` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in ``Gala`` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike. The source code for ``Gala`` has been
archived to Zenodo with the linked DOI: [@zenodo]

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge help debugging the software from Dr. Victor Francia and Mr. Kaiqiao Wu (University College London) and Kai. Also  

# References
Example paper.bib file:

@article{Pearson:2017,
    Adsnote = {Provided by the SAO/NASA Astrophysics Data System},
    Adsurl = {http://adsabs.harvard.edu/abs/2017arXiv170304627P},
    Archiveprefix = {arXiv},
    Author = {{Pearson}, S. and {Price-Whelan}, A.~M. and {Johnston}, K.~V.},
    Eprint = {1703.04627},
    Journal = {ArXiv e-prints},
    Keywords = {Astrophysics - Astrophysics of Galaxies},
    Month = mar,
    Title = {{Gaps in Globular Cluster Streams: Pal 5 and the Galactic Bar}},
    Year = 2017
}

@book{Binney:2008,
    Adsnote = {Provided by the SAO/NASA Astrophysics Data System},
    Adsurl = {http://adsabs.harvard.edu/abs/2008gady.book.....B},
    Author = {{Binney}, J. and {Tremaine}, S.},
    Booktitle = {Galactic Dynamics: Second Edition, by James Binney and Scott Tremaine.~ISBN 978-0-691-13026-2 (HB).~Published by Princeton University Press, Princeton, NJ USA, 2008.},
    Publisher = {Princeton University Press},
    Title = {{Galactic Dynamics: Second Edition}},
    Year = 2008
}

@article{zenodo,
    Abstractnote = {Gala is a Python package for Galactic astronomy and gravitational dynamics. The bulk of the package centers around implementations of gravitational potentials, numerical integration, and nonlinear dynamics.},
    Author = {Adrian Price-Whelan and Brigitta Sipocz and Syrtis Major and Semyeong Oh},
    Date-Modified = {2017-08-13 14:14:18 +0000},
    Doi = {10.5281/zenodo.833339},
    Month = {Jul},
    Publisher = {Zenodo},
    Title = {adrn/gala: v0.2.1},
    Year = {2017},
    Bdsk-Url-1 = {http://dx.doi.org/10.5281/zenodo.833339}
}

@article{gaia,
    author = {{Gaia Collaboration}},
    title = "{The Gaia mission}",
    journal = {\aap},
    archivePrefix = "arXiv",
    eprint = {1609.04153},
    primaryClass = "astro-ph.IM",
    keywords = {space vehicles: instruments, Galaxy: structure, astrometry, parallaxes, proper motions, telescopes},
    year = 2016,
    month = nov,
    volume = 595,
    doi = {10.1051/0004-6361/201629272},
    adsurl = {http://adsabs.harvard.edu/abs/2016A%26A...595A...1G},
}

@article{astropy,
    author = {{Astropy Collaboration}},
    title = "{Astropy: A community Python package for astronomy}",
    journal = {\aap},
    archivePrefix = "arXiv",
    eprint = {1307.6212},
    primaryClass = "astro-ph.IM",
    keywords = {methods: data analysis, methods: miscellaneous, virtual observatory tools},
    year = 2013,
    month = oct,
    volume = 558,
    doi = {10.1051/0004-6361/201322068},
    adsurl = {http://adsabs.harvard.edu/abs/2013A%26A...558A..33A}
}

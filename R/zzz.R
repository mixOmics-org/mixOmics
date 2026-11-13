#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created:  23-11-2015
# last modified: 18-02-2016
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname){ packageStartupMessage( "\nLoaded mixOmics ", as.character(packageDescription("mixOmics")[["Version"]]),
                                                                "\nThank you for using mixOmics!",
                                                                "\nTutorials: http://mixomics.org",
                                                                "\nBookdown vignette: https://mixomicsteam.github.io/Bookdown",
                                                                "\nQuestions, issues: Follow the prompts at http://mixomics.org/contact-us",
                                                                "\nCite us:  citation('mixOmics')\n"
                                                                
)}

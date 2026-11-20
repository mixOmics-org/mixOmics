#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created:  23-11-2015
# last modified: 18-02-2016
#
# This file is part of mixOmics.
#
# mixOmics is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, version 3.
#
# mixOmics is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
# General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with mixOmics. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: AGPL-3.0-only
#############################################################################################################

#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname){ packageStartupMessage( "\nLoaded mixOmics ", as.character(packageDescription("mixOmics")[["Version"]]),
                                                                "\nThank you for using mixOmics!",
                                                                "\nTutorials: https://mixomics.org/",
                                                                "\nVignette: https://guides.mixomics.org/mixomics-guide-vignette/",
                                                                "\nQuestions, issues: Follow the prompts at https://mixomics.org/contact/",
                                                                "\nCite us:  citation('mixOmics')\n"
                                                                
)}

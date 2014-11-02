##############################################################################
# Delta test, version 0.4.8 (beta release)
#
# Copyright (C) 2009-2013 by Maciej Eder, Jan Rybicki & Mike Kestemont.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# ############################################################################
# 
# To cite this script in publications you might use:
#   Eder, M., Rybicki, J. (2011). Stylometry with R. In "Digital Humanities 
#   2011: Conference Abstracts." Stanford University, Stanford, CA, pp. 308-11.
#
# Contact with the authors:
#     Maciej Eder <maciejeder@gmail.com>
#     Jan Rybicki <jkrybicki@gmail.com>
#     Mike Kestemont <mike.kestemont@gmail.com>
#
# ############################################################################
#
# ver. 0.4.8, 2012/12/29  --> customizable plot area, font size, etc.;
#                             thoroughly rewritten code for margins assignment;
#                             scatterplots represented either by points, or
#                             by labels, or by both (customizable label offset);
#                             saving the words (features) actually used,
#                             saving the table of actually used frequencies
# ver. 0.4.7, 2012/11/25  --> new output/input extensions: optional custom list 
#                             of files to be analyzed, saving distance table(s) 
#                             to external files, support for TXM Textometrie 
#                             Project; color cluster analysis graphs (at last!)
# ver. 0.4.6, 2012/09/09  --> code revised, cleaned, bugs fixed
# ver. 0.4.5.4, 2012/09/03 -> added 2 new PCA visualization flavors
# ver. 0.4.5.3, 2012/08/31 -> new GUI written
# ver. 0.4.5.2, 2012/08/27 -> added functionality for normal sampling
# ver. 0.4.5.1, 2012/08/22 -> support for Dutch added
# ver. 0.4.5, 2012/07/07  --> option for choosing corpus files; code cleaned;
#                             bugs fixed
# ver. 0.4.4, 2012/05/31  --> the core code rewritten, the I/II set division 
#                             abandoned, GUI remodeled, tooltips added, 
#                             different input formats supported (xml etc.), 
#                             config options loaded from external file;
#                             the code forked into (1) this script, supporting 
#                             explanatory analyses (MDS, Cons. Trees, ...),
#                             (2) a script for machine-learning methods
# ver. 0.4.3, 2012/04/28  --> feature selection (word and character n-grams)
# ver. 0.4.2, 2012/02/10  --> three ways of splitting words in English; 
#                             bugs fixed; GUI code rearranged and simplified
# ver. 0.4.1, 2011/06/27  --> better output; better text files uploading,
#                             new options for culling and ranking of candidates
# ver. 0.4.0, 2011/06/20  --> the official world-premiere (Stanford, CA)
# ver. 0.3.9b, 2011/06/1  --> the code simplified; minor cleaning
# ver. 0.3.9, 2011/05/21  --> uploading wordlist from external source;
#                             thousands of improvements; the code simplified
# ver. 0.3.8, 2010/11/01  --> skip top frequency words option added
# ver. 0.3.7, 2010/11/01  --> better graphs; attempt at better graph layout
# ver. 0.3.6, 2010/07/31  --> more graphic options; dozens of improvements
# ver. 0.3.5, 2010/07/19  --> module for color graphs; module for PCA
# ver. 0.3.4, 2010/07/12  --> module for uploading corpus files improved
# ver. 0.3.3, 2010/06/03  --> the core code simplified and improved (faster!)
# ver. 0.3.2, 2010/05/10  --> reordered GUI; minor cleaning
# ver. 0.3.1, 2010/05/10  --> the z-scores module improved
# ver. 0.3.0, 2009/12/26  --> better counter of ’good guesses’; option for
#                             randomly generated samples; minor improvements
# ver. 0.2.99, 2009/12/25 --> platform-independent outputfile saving
# ver. 0.2.98, 2009/12/24 --> GUI thoroughly integrated with initial variables
# ver. 0.2.10, 2009/11/28 --> corrected MFW display in graph, more analysis 
#                             description in outputfile
# ver. 0.2.9, 2009/11/22  --> auto graphs for MSD and CA
# ver. 0.2.8a, 2009/11/21 --> remodeled GUI
# ver. 0.2.8, 2009/11/20  --> GUI: radiobuttons, checkbuttons
# ver. 0.2.7, 2009/11/19  --> language-determined pronoun selection
# ver. 0.2.6, 2009/11/18  --> dialog box (GUI)
# ver. 0.2.5, 2009/11/16  --> module for different distance measures;
#                             thousands of improvements (I/O, interface, etc.)
# ver. 0.2.2, 2009/10/25  --> numerous little improvements; deleting pronouns
# ver. 0.2.1, 2009/08/23  --> module for culling; module for bootstrapping
# ver. 0.2.0, 2009/08/23  --> module for uploading plain text files
# ver. 0.1.9, 2009/07/31  --> numerous improvements, the code simplified
# ver. 0.1.4, 2009/07/19  --> loop for different MFW settings
# ver. 0.0.1, 2009/07/01  --> some bash and awk scripts translated into R
##############################################################################

# Before defining any new objects (variables), all the existing R objects are 
# recorded to preserve them from being removed; however, they are still 
# vulnerable to being overwritten!
variables.not.to.be.removed = ls()





label.offset = 3
add.to.margins = 2

# "labels" || "points" || "both"
text.id.on.graphs = "labels"




######################################
# COMPATIBILITY MODE options

# original algorithm for estabilishing consensus trees, used in ver. < 0.4.8;
# this is still available: if you want to replicate your old tests, say TRUE
nj.consensus.tree = FALSE

# in ver. 0.4.7, color dendrograms were available; they were produced using
# the neighbor joining (NJ) algorithm. If you want to use it, say TRUE
nj.cluster.analysis = FALSE

######################################
      # loading the required library
      if(nj.cluster.analysis == TRUE) {
        library(ape)
        }
######################################






#######  GENERAL SETTINGS (GUI/TEXT-MODE)  ###################################

# If you wish to use a simple yet effective graphical interface (GUI),
# just set the following option to TRUE, otherwise switch this option to FALSE
# and edit manually the rest of variables (see below).
# If you switch this option on, the values indicated in the following sections 
# will serve as default for the GUI for the first run of the script on a corpus. 
# In the subsequent runs, last values will appear as default in the GUI.

interactive.mode.with.GUI = TRUE


#######  TEXT- AND LANGUAGE-DEPENDENT SETTINGS  ####################

# format of corpus files; available choices are:
# "plain", "xml", "xml.drama", "xml.notitles", "html"
corpus.format = "plain"

# how many MFW ("Most frequent Words") should be taken into analysis 
# (if mfw.min value = max.mfw, then no multiple iterations will be computed)
# start.at option enables skipping top frequency words: you should
# indicate the desired start position of your list (in most cases you will 
# probably prefer setting it to 1, the rank of the single most frequent word,
# so that no words are skipped at the top of the frequency spectrum).

mfw.min = 100
mfw.max = 100
mfw.incr = 100


start.at = 1

# culling rate specifies the percentage of texts in a corpus in which a given word 
# must be found in order to be included in the analysis. Thus, a 100% culling 
# rate limits the analysis to words that appear at least once in every text 
# in the corpus; at a 50% culling rate, a word is included into the analysis 
# when it appears in at least half of the texts in the corpus; a 0% culling 
# rate (or no culling) means that no words are omitted.
# about min=max: see above

culling.min = 0
culling.max = 0
culling.incr = 20

# Deleting pronouns (this is independent of the culling procedure).
# If the "delete pronouns" option is switched to TRUE, choose one language
# of the following: English, Polish, Latin, French, German, Italian, Hungarian, Dutch, Spanish
# (the editable lists of pronouns are available below; see: advanced settings).
# Additionally, there are a few variants of language settings available:
# English.contr, English.all, and Latin.corr. Their meaning is as follows:
#     "English.contr": treats the contractions as single words, i.e. strings
#         such as "don't", "you've" etc. will not be split into two words.
#     "English.all": keeps the contractions (as above), and also prevents
#         from splitting compound words (mother-in-law, double-decker, etc.)
#     "Latin.corr": since some editions do not distinguish the letters v/u,
#         this option provides a consistent conversion to "u" in each text.

delete.pronouns = FALSE
corpus.lang = "English.all"

# Selection of features. In classical approaches, frequencies of the most
# frequent words (MFW) are used as the basis for multidimensional analyses.
# It has been argued, however, that other features are also worth considering,
# especially word and/or character n-grams. The general concept of n-grams
# is to divide a string of single words/characters into a sequence of n
# elements. Given a sample sentence "This is a simple example", the character 
# 2-grams are as follows: "th", "hi", "is", "s ", " i", "is", "s ", " a", "a ",
# " s", "si", "im", "mp", etc. The same sentence split into word 2-grams:
# "this is", "is a", "a simple", "simple sentence".
# Another question is whether it really increases the accuracy of attribution;
# further reading: Eder, M. (2011). Style-markers in authorship attribution: 
# A cross-language study of the authorial fingerprint, ’Studies in Polish
# Linguistics’ 6: 101-16.
# Two types of n-grams are available: characters (option "c"), and words ("w").

analyzed.features = "w"
ngram.size = 1

#######  MATHEMATICAL SETTINGS (DISTANCE MEASURE)  #################

# Strictly speaking, the choice of the appropriate distance measure
# is the core of the statistical procedure provided by this script.
# (However, the distance measures do not apply to the PCA method)
# Although this choice is not easy, some of the following measures
# seem to be more suitable for linguistic purposes than others.
# On theoretical grounds, Euclidean Distance and Manhattan
# Distance should be avoided in stylometry. Canberra Distance is quite 
# troublesome but effective e.g. for Latin (it should be combined with 
# careful culling settings and a limited number of MFW taken into analysis). 
# For English, usually Classic Delta is a good choice. A theoretical 
# explanation of the measures implemented in this script is pending.
#
# The available distance measures (choose ONE) are as follows:
#   "CD" --> Classic Delta as developed by Burrows
#   "AL" --> Argamon’s Linear Delta (based on Euclidean principles)
#   "ED" --> Eder’s Delta (explanation and mathematical equation: soon)
#   "ES" --> Eder’s Simple (explanation and mathematical equation: soon)
#   "MH" --> Manhattan Distance (obvious and well documented)
#   "CB" --> Canberra Distance (risky, but sometimes amazingly good)
#   "EU" --> Euclidean Distance (basic, and the most "natural")

distance.measure = "CD"


########  VISUALIZATION METHODS, LOGS, REPORTS, DISPLAY OPTIONS  ############

# Statistical method to be used; choose one:
# Cluster Analysis: "CA"
# Multidimensional Scaling: "MDS"
# Principal Components Analysis (based on a covariance table): "PCV"
# Principal Components Analysis (based on a correlation table): "PCR"
# Bootstrap Consensus Tree: "BCT".
#
# Note on the bootstrap procedure: multiple iterations will build 
# a consensus tree
#   ATTENTION: this option requires the ape library, which you can install at
# any time using the "install.packages()" command.

consensus.strength = 0.5
analysis.type="CA"


# Do you want to display the graph on the screen?
# Do you want to write the graph directly to a graphics file? Which format?
# You can display the graph on the screen AND write to a file (the latter 
# will be done with much better quality).

display.on.screen = TRUE
write.pdf.file = FALSE
write.jpg.file = FALSE
write.emf.file = FALSE    # Windows only
write.png.file = FALSE

# dimensions of the plot area (expressed in inches), font size,
# thickness of the lines used to plot the graphs.
# Since it is usually hard to remember all the values, an additional option
# is provided to reset the picture options -- if this is switched on,
# the remaining options will be overwritten

plot.options.reset = FALSE
plot.custom.height = 7
plot.custom.width = 7
plot.font.size = 10
plot.line.thickness = 1



# Do you want the graphs colored? The script will automatically assign 
# the same colors to those texts that have the same first segment of their 
# file names (the first string ending in "_").
# Available options: "colors" || "grayscale" || "black"

colors.on.graphs ="colors"

# Do you want titles on your graphs, listing the most important parameters?

titles.on.graphs = TRUE


# Layout of dendrogram: horizontal/vertical (Cluster Analysis only)

dendrogram.layout.horizontal = TRUE


# Initialize pca VISUALIZATION options; choose one:
# "classic", "loadings", "technical", "symbols"

pca.visual.flavour = "classic" # || "technical" || "symbols"


# Sometimes, you might want to save computed table(s) of distances
# for further analysis. Switch the following option TRUE to make it possible.
# The same applies to the list of words (or other features) actually used 
# in analysis, i.e. after culling, pronoun deletion, etc. Again, one might 
# want to save the table of frequencies actually used

save.distance.tables = FALSE
save.analyzed.features = FALSE
save.analyzed.freqs = FALSE


#######  ADVANCED SETTINGS (FOR EXPERTS ONLY)  ########################

# Normally, the script computes a huge table of thousands 
# of word frequencies for all texts in your corpus. This is a non-trivial 
# and time-consuming task. If done once, there is no need to waste time 
# and do it again, because the tables are also saved in the output file
# "table_with_frequencies.txt". To retrieve all the word
# frequencies from the file, switch this option to TRUE.
# BUT it MUST be set to FALSE when you switch corpora in the same R session, 
# or when you switch from word to character analysis, or change your n for
# your n-grams (or if you suddenly realize you’ve picked the
# wrong language!

use.existing.freq.tables = FALSE

# Some people like to see what’s going on, and to be able to revise/edit
# the list of words for analysis. To meet their wishes, the script
# saves the list into a separate output file, "wordlist.txt".
# You can add any words to the list and either delete as many as you want, 
# or mark the unwanted words with "#" (just like these comments are marked). 
# Switching this option on prevents the script from overwriting the file, 
# and makes sure that the wordlist is loaded from there.

use.existing.wordlist = FALSE

# Otherwise, select files manually.

interactive.files = FALSE

# Another option makes it possible to upload the files using an external list
# of files. It should be named "files_to_analyze.txt" and be put into the working
# directory. The items (i.g. file names) should be separated either by spaces,
# tabs, or newlines. The delimiters can be mixed and/or multiplied, thus even
# a very untidy list will be interpreted correctly. 

use.custom.list.of.files = TRUE

# Usually, it is recommended to cut off the tail of the word-list;
# if you do not want to cut the list, then the variable may be set to an 
# absurdly big number, or to "mfw.list.cutoff = mfw.list.of.all"
# (and then you are advised to use a fast computer).

mfw.list.cutoff = 5000

# pronouns (and other words) to be deleted
# * what are the selection criteria used here? Personal, possessive, ...? *

# Polish
pol.pronouns = c("ci", "ciebie", "cię", "go", "ich", "im", "ja", "ją", "je", "jego", "jej", "jemu", "ma", "mą", "me", "mego", "mej", "memu", "mi", "mną", "mnie", "moi", "moich", "moim", "moimi", "moja", "moją", "moje", "mojego", "mojej", "mojemu", "mój", "mu", "my", "mych", "mym", "mymi", "nam", "nami", "nas", "nią", "nich", "nie", "niego", "niej", "niemu", "nim", "nimi", "on", "ona", "one", "oni", "ono", "swa", "swą", "swe", "swego", "swej", "swemu", "swoi", "swoich", "swoim", "swoimi", "swoja", "swoją", "swoje", "swojego", "swojej", "swojemu", "swój", "swych", "swym", "swymi", "tobą", "tobie", "twa", "twą", "twe", "twego", "twej", "twemu", "twoi", "twoich", "twoim", "twoimi", "twoja", "twoją", "twoje", "twojego", "twojej", "twojemu", "twój", "twych", "twym", "twymi", "ty", "wam", "wami", "was", "wy", "wasz", "wasza", "wasze", "waszym", "waszymi", "waszych", "waszego", "waszej", "waszą")
# English
eng.pronouns = c("he", "her", "hers", "herself", "him", "himself", "his", "i", "me", "mine", "my", "myself", "our", "ours", "ourselves", "she", "thee", "their", "them", "themselves", "they", "thou", "thy", "thyself", "us", "we", "ye", "you", "your", "yours", "yourself")
# Latin
lat.pronouns = c("ea", "eae", "eam", "earum", "eas", "ego", "ei", "eis", "eius", "eo", "eorum", "eos", "eum", "id", "illa", "illae", "illam", "illarum", "illas", "ille", "illi", "illis", "illius", "illo", "illorum", "illos", "illud", "illum", "is", "me", "mea", "meae", "meam", "mearum", "meas", "mei", "meis", "meo", "meos", "meorum", "meum", "meus", "mihi", "nobis", "nos", "noster", "nostra", "nostrae", "nostram", "nostrarum", "nostras", "nostri", "nostris", "nostro", "nostros", "nostrorum", "nostrum", "sua", "suae", "suam", "suarum", "suas", "sui", "suis", "suo", "suos", "suorum", "suum", "suus", "te", "tibi", "tu", "tua", "tuae", "tuam", "tuarum", "tuas", "tui", "tuis", "tuo", "tuos", "tuorum", "tuum", "tuus", "vester", "vestra", "vestrae", "vestram", "vestrarum", "vestras", "vestri", "vestris", "vestro", "vestros", "vestrorum", "vestrum", "vobis", "vos")
# French
fra.pronouns = c("je", "me", "moi", "tu", "te", "toi", "il", "elle", "le", "la", "lui", "se", "lui", "elle", "soi", "nous", "vous", "ils", "elles", "les", "leur", "se", "eux", "elles", "soi")
# German
ger.pronouns = c("ich", "mich", "mir", "mein", "meine", "meiner", "meines", "du", "dich", "dir", "dein", "deine", "deiner", "deines", "er", "sich", "ihr", "ihrer", "ihn", "ihnen", "sein", "seiner", "seines", "seine", "sie", "wir", "uns", "unser", "unsere", "euch", "eure", "euer")
# Italian
ita.pronouns = c("ci", "gli", "io", "la", "le", "lei", "li", "loro", "lo", "lui", "me", "mi", "noi", "si", "te", "ti", "tu", "vi", "voi", "egli", "ella", "esso", "essa", "essi", "esse", "mio", "mia", "miei", "mie", "tuo", "tua", "tuoi", "tue", "suo", "sua", "suoi", "sue", "nostro", "nostra", "nostri", "nostre", "vostro", "vostra", "vostri", "vostre", "loro", "loro", "loro", "loro")
# Hungarian
hun.pronouns = c("annak", "az", "azzal", "bele", "belé", "beléd", "beléje", "beléjük", "belém", "belénk", "belétek", "belöle", "belőled", "belőlem", "belőletek", "belőlük", "belőlünk", "benne", "benned", "bennem", "bennetek", "bennük", "bennünk", "én", "ennek", "enyéim", "enyém", "enyémek", "érte", "érted", "értem", "értetek", "értük", "értünk", "ez", "ezzel", "hozzá", "hozzád", "hozzája", "hozzájuk", "hozzám", "hozzánk", "hozzátok", "maga", "magáé", "magáéi", "maguk", "maguké", "magukéi", "mi", "mieink", "mienk", "miénk", "nála", "nálad", "nálam", "nálatok", "náluk", "nálunk", "neked", "nekem", "neki", "nekik", "nektek", "nekünk", "ő", "ők", "ön", "öné", "önéi", "önnek", "önnel", "önök", "önöké", "önökéi", "önökkel", "önöknek", "övé", "övéi", "övéik", "övék", "rád", "rája", "rajta", "rajtad", "rajtam", "rajtatok", "rajtuk", "rajtunk", "rájuk", "rám", "ránk", "rátok", "róla", "rólad", "rólam", "rólatok", "róluk", "rólunk", "te", "ti", "tied", "tiéd", "tieid", "tieitek ", "tietek", "tiétek", "tőle", "tőled", "tőlem", "töletek", "tőlük", "tőlünk", "vele", "veled", "velem", "veletek", "velük", "velünk")
# Dutch
dut.pronouns = c("hij", "haar", "haarzelf", "hijzelf", "hemzelf", "hem", "ik", "ikzelf", "mijn", "mij", "mijzelf", "me", "mezelf", "zich", "zichzelf", "ons", "onze", "onszelf", "u", "uw", "uzelf", "zij", "zijzelf", "wij", "wijzelf", "jij", "jijzelf", "jouw", "jouwe", "jou", "jouzelf", "elkaar", "hen", "henzelf", "hun", "hunzelf", "zich", "elkaar", "wie", "wat", "welke")
# Spanish
sp.pronouns = c("yo", "me", "mí", "tú", "te", "ti", "usted", "ud", "le", "lo", "la", "se", "sí", "él", "lo", "ella", "nos", "nosotros", "nosotras", "vosotros", "vosotras", "ustedes", "ud", "les", "los", "las", "se", "ellos", "los", "ellas")



# This option enables integration with TXM corpus management system
# (see: Textometrie Project, http://textometrie.ens-lyon.fr/).
# Normally, it should be switched off, since it is used only when the script
# is invoked from inside the TXM environment. WARNING: experimental.

txm.compatibility.mode = FALSE



########  SAMPLING OPTIONS  ############


sampling = "no.sampling" # || "random.sampling" || "no.sampling"

# When dealing with longer text, one might want to divide these in samples of 
# an equal size. This can be achieved by setting the sampling variable
# (default="normal.sampling") and specifying the number of words per sample 
# via the sample.size parameter: Integer, default=10000).

sample.size = 10000 # expressed in words, also if you’re using character n-grams

# when the analyzed texts are significantly unequal in length, it is not a bad
# idea to prepare samples as randomly chosen "bags of words". For this, set the
# sampling variable to "random.sampling". The desired size of the sample should
# be indicated via the length.of.random.sample variable.
# Sampling with and without replacement is also available.
# (Further reading: Eder, M. (2010). Does Size Matter? Authorship Attribution,
# Short Samples, Big Problem. In "Digital Humanities 2010: Conference 
# Abstracts." King’s College London 2010, pp. 132-35.)
#
# ATTENTION: this makes sense only if "use.existing.freq.tables" is set "FALSE"

length.of.random.sample = 10000
sampling.with.replacement = FALSE

# It is also possible to use the entire corpus texts as samples (regardless 
# of their length and differences therein). For this, set the sampling variable 
# to "no.sampling"


# the variables are now ready to use (unless the GUI option was chosen)
# ###################################################################



# #################################################
# sanity check for some of the initial variables -- just in case
# #################################################

# Given a language option ("English", "Polish", "Latin" etc., as described 
# above), this procedure selects one of the lists of pronouns
# If no language was chosen (or if a desired language is not supported, or if 
# there was a spelling mistake), then the variable will be set to "English". 
# If "Pronouns deleted" is set to FALSE, this is immaterial.

if(exists("pronouns") == FALSE) # checking if the "pronouns" box is empty
    pronouns = eng.pronouns

# This prevents us from choosing a non-existing distance measure -- in such
# case the default distance (Classic Delta) will be switched on. Be aware
# of correct spelling: then the default value will be assigned as well!

if(distance.measure %in% c("CD","AL","ED","ES","MH","CB","EU") == FALSE) {
  distance.measure = "CD"
}


# #################################################



# #################################################
#
# the GUI module 
#
# #################################################

# At the beginning of the script, you could decide whether use the GUI module 
# or not; if the appropriate option was switched on, the GUI will start now;
# Since it’s written in TclTk, with some additional twists, you need to install
# the tcltk2 package (on top of the regular tcltk, which is usually installed 
# with R anyway.

if (interactive.mode.with.GUI == TRUE) {
  library(tcltk)
  library(tcltk2)

if(file.exists("config.txt") == TRUE) {
  source("config.txt") }

# ###################################################################

.Tcl("font create myDefaultFont -family tahoma -size 8")
.Tcl("option add *font myDefaultFont")  

  cancel_pause <- FALSE
  tt <- tktoplevel()
  tktitle(tt) <- "Stylometry with R: enter analysis parameters"
  
  push_OK <- function(){
      cancel_pause <<- TRUE
      tkdestroy(tt)
      }

corpus.format <- tclVar(corpus.format)
mfw.min <- tclVar(mfw.min)
mfw.max <- tclVar(mfw.max)
mfw.incr <- tclVar(mfw.incr)
start.at <- tclVar(start.at)
culling.min <- tclVar(culling.min)
culling.max <- tclVar(culling.max)
culling.incr <- tclVar(culling.incr)
ngram.size <- tclVar(ngram.size)
analyzed.features <- tclVar(analyzed.features)
use.existing.freq.tables <- tclVar(use.existing.freq.tables)
use.existing.wordlist <- tclVar(use.existing.wordlist)
interactive.files <- tclVar(interactive.files)
use.custom.list.of.files <- tclVar(use.custom.list.of.files)
mfw.list.cutoff <- tclVar(mfw.list.cutoff)
analysis.type <- tclVar(analysis.type)
delete.pronouns <- tclVar(delete.pronouns)
corpus.lang <- tclVar(corpus.lang)
distance.measure <- tclVar(distance.measure)
display.on.screen <- tclVar(display.on.screen)
write.pdf.file <- tclVar(write.pdf.file)
write.jpg.file <- tclVar(write.jpg.file)
write.emf.file <- tclVar(write.emf.file)
write.png.file <- tclVar(write.png.file)
colors.on.graphs <- tclVar(colors.on.graphs)
titles.on.graphs <- tclVar(titles.on.graphs)
dendrogram.layout.horizontal <- tclVar(dendrogram.layout.horizontal)
pca.visual.flavour <- tclVar(pca.visual.flavour)
save.distance.tables <- tclVar(save.distance.tables)
save.analyzed.features <- tclVar(save.analyzed.features)
save.analyzed.freqs <- tclVar(save.analyzed.freqs)
sampling <- tclVar(sampling)
sample.size <- tclVar(sample.size)
length.of.random.sample <- tclVar(length.of.random.sample)
consensus.strength <- tclVar(consensus.strength)
plot.options.reset <- tclVar(plot.options.reset)
plot.custom.height <- tclVar(plot.custom.height)
plot.custom.width <- tclVar(plot.custom.width)
plot.font.size <- tclVar(plot.font.size)
plot.line.thickness <- tclVar(plot.line.thickness)
text.id.on.graphs <- tclVar(text.id.on.graphs)
add.to.margins <- tclVar(add.to.margins)
label.offset <- tclVar(label.offset)






f1 <- tkframe(tt)
f2 <- tkframe(tt)
f3 <- tkframe(tt)
f4 <- tkframe(tt)
f5 <- tkframe(tt)

# layout of the GUI begins here:
tab1 <- function() {
tkgrid(f1,row=1,column=0,columnspan=5)
tkgrid.forget(f2)
tkgrid.forget(f3)
tkgrid.forget(f4)
tkgrid.forget(f5)
tkconfigure(t1.but,state="disabled", background="white")
tkconfigure(t2.but,state="normal", background="aliceblue")
tkconfigure(t3.but,state="normal", background="aliceblue")
tkconfigure(t4.but,state="normal", background="aliceblue")
tkconfigure(t5.but,state="normal", background="aliceblue")
}
tab2 <- function() {
tkgrid(f2,row=1,column=0,columnspan=5)
tkgrid.forget(f1)
tkgrid.forget(f3)
tkgrid.forget(f4)
tkgrid.forget(f5)
tkconfigure(t2.but,state="disabled", background="white")
tkconfigure(t1.but,state="normal", background="aliceblue")
tkconfigure(t3.but,state="normal", background="aliceblue")
tkconfigure(t4.but,state="normal", background="aliceblue")
tkconfigure(t5.but,state="normal", background="aliceblue")
}
tab3 <- function() {
tkgrid(f3,row=1,column=0,columnspan=5)
tkgrid.forget(f1)
tkgrid.forget(f2)
tkgrid.forget(f4)
tkgrid.forget(f5)
tkconfigure(t3.but,state="disabled", background="white")
tkconfigure(t1.but,state="normal", background="aliceblue")
tkconfigure(t2.but,state="normal", background="aliceblue")
tkconfigure(t4.but,state="normal", background="aliceblue")
tkconfigure(t5.but,state="normal", background="aliceblue")
}
tab4 <- function() {
tkgrid(f4,row=1,column=0,columnspan=5)
tkgrid.forget(f1)
tkgrid.forget(f2)
tkgrid.forget(f3)
tkgrid.forget(f5)
tkconfigure(t4.but,state="disabled", background="white")
tkconfigure(t1.but,state="normal", background="aliceblue")
tkconfigure(t2.but,state="normal", background="aliceblue")
tkconfigure(t3.but,state="normal", background="aliceblue")
tkconfigure(t5.but,state="normal", background="aliceblue")
}
tab5 <- function() {
tkgrid(f5,row=1,column=0,columnspan=5)
tkgrid.forget(f1)
tkgrid.forget(f2)
tkgrid.forget(f3)
tkgrid.forget(f4)
tkconfigure(t5.but,state="disabled", background="white")
tkconfigure(t1.but,state="normal", background="aliceblue")
tkconfigure(t2.but,state="normal", background="aliceblue")
tkconfigure(t3.but,state="normal", background="aliceblue")
tkconfigure(t4.but,state="normal", background="aliceblue")
}
t1.but <- tkbutton(tt,text="     INPUT & LANGUAGE     ",command=tab1)
t2.but <- tkbutton(tt,text="         FEATURES         ",command=tab2)
t3.but <- tkbutton(tt,text="        STATISTICS        ",command=tab3)
t4.but <- tkbutton(tt,text="         SAMPLING         ",command=tab4)
t5.but <- tkbutton(tt,text="          OUTPUT          ",command=tab5)
tkgrid(t1.but)
tkgrid(t2.but, column=1, row=0)
tkgrid(t3.but, column=2, row=0)
tkgrid(t4.but, column=3, row=0)
tkgrid(t5.but, column=4, row=0)
# Grid for individual tabs
 
# initial state!
tkgrid(f1,row=1,column=0,columnspan=5)
tkconfigure(t1.but,state="disabled", background="white")
tkconfigure(t2.but,state="normal", background="aliceblue")
tkconfigure(t3.but,state="normal", background="aliceblue")
tkconfigure(t4.but,state="normal", background="aliceblue")
tkconfigure(t5.but,state="normal", background="aliceblue")

# the OK button: active on each tab
#
button_1 <- tkbutton(tt,text="       OK       ",command=push_OK,relief="raised",background="aliceblue")
tkbind(button_1,"<Return>",push_OK) 
tkgrid(button_1,columnspan=10)
tk2tip(button_1, "Press this only if you've visited all the tabs, or if you know\nyou want to leave values in some as they are.")


########################################################################################################################
# layout of the GUI begins here:
#
tkgrid(tklabel(f1,text="    "),padx=0,pady=0) # blank line (serving as the top margin)
tkgrid(tklabel(f2,text="    ")) # blank line (serving as the top margin)
tkgrid(tklabel(f3,text="    ")) # blank line (serving as the top margin)
tkgrid(tklabel(f4,text="    ")) # blank line (serving as the top margin)
tkgrid(tklabel(f5,text="    ")) # blank line (serving as the top margin)

# first row: INPUT
#
entry_TXT <- tkradiobutton(f1)
entry_XML <- tkradiobutton(f1)
entry_XMLDrama <- tkradiobutton(f1)
entry_XMLNoTitles <- tkradiobutton(f1)
entry_HTML <- tkradiobutton(f1)
#
tkconfigure(entry_TXT,variable=corpus.format,value="plain")
tkconfigure(entry_XML,variable=corpus.format,value="xml")
tkconfigure(entry_XMLDrama,variable=corpus.format,value="xml.drama")
tkconfigure(entry_XMLNoTitles,variable=corpus.format,value="xml.notitles")
tkconfigure(entry_HTML,variable=corpus.format,value="html")
#
entrylabel_TXT <- tklabel(f1,text="plain text")
entrylabel_XML <- tklabel(f1,text="xml")
entrylabel_XMLDrama <- tklabel(f1,text="xml (plays)")
entrylabel_XMLNoTitles <- tklabel(f1,text="xml (no titles)")
entrylabel_HTML <- tklabel(f1,text="html")
#
tkgrid(tklabel(f1,text="       INPUT:"),entrylabel_TXT,entrylabel_XML,entrylabel_XMLDrama,entrylabel_XMLNoTitles,entrylabel_HTML,columnspan=1)
tkgrid(tklabel(f1,text="            "),entry_TXT,entry_XML,entry_XMLDrama,entry_XMLNoTitles,entry_HTML,columnspan=1)
# Tooltips for the above
tk2tip(entrylabel_TXT, "Plain text files. \nIf your corpus does not contain diacritics, no encoding is needed. \nOtherwise, use ANSI for Windows, UTF-8 for Mac/Linux.")
tk2tip(entrylabel_XML, "XML: all tags and TEI headers are removed.")
tk2tip(entrylabel_XMLDrama, "XML for plays: all tags, TEI headers, \nand speakers' names between <speaker>...</speaker> tags are removed.")
tk2tip(entrylabel_XMLNoTitles, "XML contents only: all tags, TEI headers, \nand chapter/section (sub)titles between <head>...</head> tags are removed.")
tk2tip(entrylabel_HTML, "HTML headers, menus, links and other tags are removed.")
tkgrid(tklabel(f1,text="    ")) # blank line for aesthetic purposes

# next row: LANGUAGE
#
entry_ENG <- tkradiobutton(f1)
entry_EN2 <- tkradiobutton(f1)
entry_EN3 <- tkradiobutton(f1)
entry_POL <- tkradiobutton(f1)
entry_LAT <- tkradiobutton(f1)
entry_LA2 <- tkradiobutton(f1)
entry_FRA <- tkradiobutton(f1)
entry_GER <- tkradiobutton(f1)
entry_HUN <- tkradiobutton(f1)
entry_ITA <- tkradiobutton(f1)
entry_DUT <- tkradiobutton(f1)
entry_SPA <- tkradiobutton(f1)
#
tkconfigure(entry_ENG,variable=corpus.lang,value="English")
tkconfigure(entry_EN2,variable=corpus.lang,value="English.contr")
tkconfigure(entry_EN3,variable=corpus.lang,value="English.all")
tkconfigure(entry_LAT,variable=corpus.lang,value="Latin")
tkconfigure(entry_LA2,variable=corpus.lang,value="Latin.corr")
tkconfigure(entry_POL,variable=corpus.lang,value="Polish")
tkconfigure(entry_FRA,variable=corpus.lang,value="French")
tkconfigure(entry_GER,variable=corpus.lang,value="German")
tkconfigure(entry_HUN,variable=corpus.lang,value="Hungarian")
tkconfigure(entry_ITA,variable=corpus.lang,value="Italian")
tkconfigure(entry_DUT,variable=corpus.lang,value="Dutch")
tkconfigure(entry_SPA,variable=corpus.lang,value="Spanish")
#
entrylabel_ENG <- tklabel(f1,text="    English     ")
entrylabel_POL <- tklabel(f1,text="    Polish      ")
entrylabel_LAT <- tklabel(f1,text="    Latin       ")
entrylabel_FRA <- tklabel(f1,text="    French      ")
entrylabel_GER <- tklabel(f1,text="    German      ")
entrylabel_HUN <- tklabel(f1,text="   Hungarian    ")
entrylabel_ITA <- tklabel(f1,text="    Italian     ")
entrylabel_EN2 <- tklabel(f1,text="English (contr.)")
entrylabel_EN3 <- tklabel(f1,text="  English (ALL) ")
entrylabel_LA2 <- tklabel(f1,text="Latin (u/v > u) ")
entrylabel_DUT <- tklabel(f1,text="     Dutch      ")
entrylabel_SPA <- tklabel(f1,text="    Spanish     ")
#
tkgrid(tklabel(f1,text="LANGUAGE: "),entrylabel_ENG,entrylabel_EN2,entrylabel_EN3,entrylabel_LAT,entrylabel_LA2)
tkgrid(tklabel(f1,text="          "),entry_ENG,entry_EN2,entry_EN3,entry_LAT,entry_LA2)
tkgrid(tklabel(f1,text="          "),entrylabel_POL,entrylabel_HUN,entrylabel_FRA,entrylabel_ITA,entrylabel_SPA)
tkgrid(tklabel(f1,text="          "),entry_POL,entry_HUN,entry_FRA,entry_SPA,entry_ITA)
tkgrid(tklabel(f1,text="          "),entrylabel_DUT,entrylabel_GER)
tkgrid(tklabel(f1,text="          "),entry_DUT,entry_GER)
tkgrid(tklabel(f1,text="    ")) # blank line for aesthetic purposes

# Tooltips for the above
tk2tip(entrylabel_ENG, "Plain English: contractions and \ncompound words are split")
tk2tip(entrylabel_POL, "Plain Polish: contractions and \ncompound words are split")
tk2tip(entrylabel_EN2, "Modified English: \ncontractions are not split")
tk2tip(entrylabel_EN3, "Further Modified English: contractions \nand compound words are not split")
tk2tip(entrylabel_LAT, "Plain Latin: U and V \ntreated as distinct letters")
tk2tip(entrylabel_FRA, "Plain French: contractions and \ncompound words are split")
tk2tip(entrylabel_GER, "Plain German: contractions and \ncompound words are split")
tk2tip(entrylabel_HUN, "Plain Hungarian: contractions and \ncompound words are split")
tk2tip(entrylabel_ITA, "Plain Italian: contractions and \ncompound words are split")
tk2tip(entrylabel_LA2, "Modified Latin: U and V \nboth treated as U")
tk2tip(entrylabel_DUT, "Plain Dutch: contractions and \ncompound words are split")
tk2tip(entrylabel_SPA, "Plain Castilian: contractions and \ncompound words are split")

# next row: TEXT FEATURES
entry_W <- tkradiobutton(f2)
entry_L <- tkradiobutton(f2)
cb_NGRAMS <- tkcheckbutton(f2)
entry_NGRAMSIZE <- tkentry(f2,textvariable=ngram.size,width="8")
#
tkconfigure(entry_W,variable=analyzed.features,value="w")
tkconfigure(entry_L,variable=analyzed.features,value="c")
#
entrylabel_W <- tklabel(f2,text="words")
entrylabel_L <- tklabel(f2,text="chars")
entrylabel_NGRAMSIZE <- tklabel(f2,text="ngram size")
#
tkgrid(tklabel(f2,text="        FEATURES:"),entrylabel_W,entrylabel_L,entrylabel_NGRAMSIZE)
tkgrid(tklabel(f2,text="                 "),entry_W,entry_L,entry_NGRAMSIZE)

# Tooltips for the above
tk2tip(entrylabel_W, "Select this to work on words")
tk2tip(entrylabel_L, "Select this to work on characters \n(does not make much sense unless you use ngrams)")
tk2tip(entrylabel_NGRAMSIZE, "State your n for n-grams \nto work on word/char clusters of n")
tkgrid(tklabel(f2,text="    ")) # blank line for aesthetic purposes

# next row: MFW SETTINGS
#
entry_MFW_MIN <- tkentry(f2,textvariable=mfw.min,width="8")
entry_MFW_MAX <- tkentry(f2,textvariable=mfw.max,width="8")
entry_MFW_INCR <- tkentry(f2,textvariable=mfw.incr,width="8")
entry_START_AT <- tkentry(f2,textvariable=start.at,width="8")
#
entrylabel_MFW_MIN <- tklabel(f2,text="Minimum")
entrylabel_MFW_MAX <- tklabel(f2,text="Maximum")
entrylabel_MFW_INCR <- tklabel(f2,text="Increment")
entrylabel_START_AT <- tklabel(f2,text="Start at freq. rank")
#
tkgrid(tklabel(f2,text="MFW SETTINGS:"),entrylabel_MFW_MIN,entrylabel_MFW_MAX,entrylabel_MFW_INCR,entrylabel_START_AT)
tkgrid(tklabel(f2,text="             "),entry_MFW_MIN,entry_MFW_MAX,entry_MFW_INCR,entry_START_AT)
tkgrid(tklabel(f2,text="    ")) # blank line for aesthetic purposes

# Tooltips for the above
tk2tip(entrylabel_MFW_MIN, "Set the minimum number of most frequent words. \nThe script will conduct its first analysis for \nthe number of words specified here")
tk2tip(entrylabel_MFW_MAX, "Set the maximum number of most frequent words. \nThe script will conduct its final analysis for \nthe number of words specified here")
tk2tip(entrylabel_MFW_INCR, "Set the increment added to \nthe minimum number of most frequent \nwords for each subsequent analysis.")
tk2tip(entrylabel_START_AT, "Set the number of words from the top of \nthe frequency list to skip in the analysis.")

# next row: CULLING
#
cb_DEL_PRON <- tkcheckbutton(f2)
#
entry_CUL_MIN <- tkentry(f2,textvariable=culling.min,width="8")
entry_CUL_MAX <- tkentry(f2,textvariable=culling.max,width="8")
entry_CUL_INCR <- tkentry(f2,textvariable=culling.incr,width="8")
entry_CUT_OFF <- tkentry(f2,textvariable=mfw.list.cutoff,width="8")
tkconfigure(cb_DEL_PRON,variable=delete.pronouns)
#
entrylabel_CUL_MIN <- tklabel(f2,text="Minimum")
entrylabel_CUL_MAX <- tklabel(f2,text="Maximum")
entrylabel_CUL_INCR <- tklabel(f2,text="Increment")
entrylabel_CUT_OFF <- tklabel(f2,text="List Cutoff")
cblabel_DEL_PRON <- tklabel(f2,text="Delete pronouns")
#
tkgrid(tklabel(f2,text="         CULLING:"),entrylabel_CUL_MIN,entrylabel_CUL_MAX, entrylabel_CUL_INCR,entrylabel_CUT_OFF,cblabel_DEL_PRON)
tkgrid(tklabel(f2,text="                 "),entry_CUL_MIN,entry_CUL_MAX,entry_CUL_INCR,entry_CUT_OFF,cb_DEL_PRON)
tkgrid(tklabel(f2,text="    ")) # blank line for aesthetic purposes
  
# next row: LISTS & FILES
#
cb_FREQS <- tkcheckbutton(f2)
cb_LISTS <- tkcheckbutton(f2)
cb_INTFILES <- tkcheckbutton(f2)
cb_MYFILES <- tkcheckbutton(f2)
#
tkconfigure(cb_FREQS,variable=use.existing.freq.tables)
tkconfigure(cb_LISTS,variable=use.existing.wordlist)
tkconfigure(cb_INTFILES,variable=interactive.files)
tkconfigure(cb_MYFILES,variable=use.custom.list.of.files)
#
cblabel_FREQS <- tklabel(f2,text="Existing frequencies")
cblabel_LISTS <- tklabel(f2,text="Existing wordlist")
cblabel_INTFILES <- tklabel(f2,text="Select files manually")
cblabel_MYFILES <- tklabel(f2,text="List of files")
#
tkgrid(tklabel(f2,text="    VARIOUS:"),cblabel_FREQS,cblabel_LISTS,cblabel_INTFILES,cblabel_MYFILES)
tkgrid(tklabel(f2,text="            "),cb_FREQS,cb_LISTS,cb_INTFILES,cb_MYFILES)
  
# Tooltips for the above  
tk2tip(entrylabel_CUL_MIN, "State the minimum culling setting. \n0 means no words are omitted from the analysis. \n50 means a word needs to appear in \nat least 50% of the texts to be included in the analysis. \n100 means that only words appearing in all the texts \nwill be included in the analysis")
tk2tip(entrylabel_CUL_MAX, "State the maximum culling setting. \n0 means no words are omitted from the analysis. \n50 means a word needs to appear in \nat least 50% of the texts to be included in the analysis. \n100 means that only words appearing in all the texts \nwill be included in the analysis")
tk2tip(entrylabel_CUL_INCR, "State the increment added to the minimum culling \nsetting for each subsequent analysis.")
tk2tip(entrylabel_CUT_OFF, "Set the maximum size of the word frequency table. \nAnything above 5000 requires patience and a fast computer")
tk2tip(cblabel_DEL_PRON, "Select if you want to omit pronouns in the analysis. \nThis improves attribution in some languages")
tk2tip(cblabel_FREQS, "Select to use the frequency lists generated by the previous analysis. \nThis speeds up the process dramatically. \nA very bad idea if you've just changed your selection of texts!")
tk2tip(cblabel_LISTS, "Select to use the wordlist generated by \nthe previous analysis or a custom wordlist.")
tk2tip(cblabel_INTFILES, "Select this to manually select files \nrather than use the entire corpus. \nMake sure that \"Existing frequencies\" is unchecked!")
tk2tip(cblabel_MYFILES, "Select this if you want to use custom list \nof files to be loaded; the list should be stored \nin \"files_to_analyze.txt\", and the entries delimited \nwith spaces, tabs, or newlines.")


tkgrid(tklabel(f2,text="    ")) # blank line for aesthetic purposes 

# next row: STATISTICS
#
entry_CA <- tkradiobutton(f3)
entry_MDS <- tkradiobutton(f3)
entry_PCA1 <- tkradiobutton(f3)
entry_PCA2 <- tkradiobutton(f3)
entry_CONS_TREE <- tkradiobutton(f3)
entry_CONSS <- tkentry(f3,textvariable=consensus.strength,width="8")
#
tkconfigure(entry_CA,variable=analysis.type,value="CA") # cluster.analysis
tkconfigure(entry_MDS,variable=analysis.type,value="MDS") # multidimensional.scaling
tkconfigure(entry_PCA1,variable=analysis.type,value="PCV") # pca.covariance.table
tkconfigure(entry_PCA2,variable=analysis.type,value="PCR") # pca.correlation.table
tkconfigure(entry_CONS_TREE,variable=analysis.type,value="BCT") # make.consensus.tree

#
entrylabel_CA <- tklabel(f3,text="Cluster Analysis")
entrylabel_MDS <- tklabel(f3,text="MDS")
entrylabel_PCA1 <- tklabel(f3,text="PCA (cov.)")
entrylabel_PCA2 <- tklabel(f3,text="PCA (corr.)")
entrylabel_CONS_TREE <- tklabel(f3,text="Consensus Tree")
entrylabel_CONSS <- tklabel(f3,text="Consensus strength")
#
tkgrid(tklabel(f3,text=" STATISTICS:"),entrylabel_CA,entrylabel_MDS,entrylabel_PCA1,entrylabel_PCA2)
tkgrid(tklabel(f3,text="            "),entry_CA,entry_MDS,entry_PCA1,entry_PCA2)
tkgrid(tklabel(f3,text="            "),entrylabel_CONS_TREE,entrylabel_CONSS)
tkgrid(tklabel(f3,text="            "),entry_CONS_TREE,entry_CONSS)

# Tooltips for the above
tk2tip(entrylabel_CA, "Select to perform Cluster Analysis of Delta distance table. \nThis only makes sense if there is a single iteration \n(or only a few), so set MFW_MIN and MFW_MAX \nto equal values, which in turn makes the MFW_INCR setting immaterial. \nThen do the same for your culling settings.")
tk2tip(entrylabel_MDS, "Select to perform Multidimensional Scaling of Delta distance table. \nThis only makes sense if there is a single iteration \n(or only a few), so set MFW_MIN and MFW_MAX \nto equal values, which in turn makes the MFW_INCR setting immaterial. \nThen do the same for your culling settings.")
tk2tip(entrylabel_PCA1, "Select to perform Principal Components Analysis based on a covariance matrix of Delta distance table. \nThis only makes sense if there is a single iteration (or only a few), so set MFW_MIN and MFW_MAX \nto equal values, which in turn makes the MFW_INCR setting immaterial. \nThen do the same for your culling settings.")
tk2tip(entrylabel_PCA2, "Select to perform Principal Components Analysis based on a correlation matrix of Delta distance table. \nThis only makes sense if there is a single iteration (or only a few), so set MFW_MIN and MFW_MAX \nto equal values, which in turn makes the MFW_INCR setting immaterial. \nThen do the same for your culling settings.")
tk2tip(entrylabel_CONS_TREE, "Select to perform multiple iterations of Cluster Analysis of Delta distance table \nresulting in a Bootstrap COnsensus Tree. This only makes sense \nif you have at least three valid iterations, so set MFW_MIN and MFW_MAX, \nand/or CUL_MIN and CUL_MAX to different values.")
tk2tip(entrylabel_CONSS, "Select to set the consensus strength for the Bootstrap Tree. \nOnly makes sense if you select that option is checked above. \nLegal values are from 0.4 (40% underlying CA graphs need to agree \non a given connection) to 1 (all underlying CA graphs need to agree).")

# next row: DISTANCES
#
entry_CD <- tkradiobutton(f3)
entry_AL <- tkradiobutton(f3)
entry_ED <- tkradiobutton(f3)
entry_ES <- tkradiobutton(f3)
entry_MH <- tkradiobutton(f3)
entry_CB <- tkradiobutton(f3)

entry_EU <- tkradiobutton(f3)
#
tkconfigure(entry_CD,variable=distance.measure,value="CD")
tkconfigure(entry_AL,variable=distance.measure,value="AL")
tkconfigure(entry_ED,variable=distance.measure,value="ED")
tkconfigure(entry_ES,variable=distance.measure,value="ES")
tkconfigure(entry_MH,variable=distance.measure,value="MH")
tkconfigure(entry_CB,variable=distance.measure,value="CB")
tkconfigure(entry_EU,variable=distance.measure,value="EU")
#
entrylabel_CD <- tklabel(f3,text="Classic Delta")
entrylabel_AL <- tklabel(f3,text="Argamon's Delta")
entrylabel_ED <- tklabel(f3,text="Eder's Delta")
entrylabel_ES <- tklabel(f3,text="Eder's Simple")
entrylabel_MH <- tklabel(f3,text="Manhattan")
entrylabel_CB <- tklabel(f3,text="Canberra")
entrylabel_EU <- tklabel(f3,text="Euclidean")
#
tkgrid(tklabel(f3,text="  DISTANCES:"),entrylabel_CD,entrylabel_AL,entrylabel_ED,entrylabel_ES)
tkgrid(tklabel(f3,text="            "),entry_CD,entry_AL,entry_ED,entry_ES)
tkgrid(tklabel(f3,text="            "),entrylabel_MH,entrylabel_CB,entrylabel_EU)
tkgrid(tklabel(f3,text="            "),entry_MH,entry_CB,entry_EU)
tkgrid(tklabel(f3,text="    ")) # blank line for aesthetic purposes

# Tooltips for the above
tk2tip(entrylabel_CD, "Select the Classic Delta measure as developed by Burrows.")
tk2tip(entrylabel_AL, "Select Argamon's Linear Delta (based on Euclidean principles).")
tk2tip(entrylabel_ED, "Select Eder's Delta (explanation and mathematical equation: TBA).")
tk2tip(entrylabel_ES, "Select Eder's Simple measure (explanation and mathematical equation: TBA).")
tk2tip(entrylabel_MH, "Select Manhattan Distance (obvious and well documented).")
tk2tip(entrylabel_CB, "Select Canberra Distance (risky, but sometimes amazingly good).")
tk2tip(entrylabel_EU, "Select Euclidean Distance (basic and the most *natural*).")

# next row: SAMPLING
entry_SAMP <- tkradiobutton(f4)
entry_RAND <- tkradiobutton(f4)
entry_NOSAMP <- tkradiobutton(f4)
  
tkconfigure(entry_SAMP, variable=sampling, value="normal.sampling")
tkconfigure(entry_RAND, variable=sampling, value="random.sampling")
tkconfigure(entry_NOSAMP, variable=sampling, value="no.sampling")
  
entry_SAMPLESIZE <- tkentry(f4,textvariable=sample.size,width="10")
entry_SIZE <- tkentry(f4,textvariable=length.of.random.sample,width="10")

entrylabel_SAMP <- tklabel(f4,text="Normal sampling")
entrylabel_RAND <- tklabel(f4,text="Random sampling")
entrylabel_NOSAMP <- tklabel(f4,text="No sampling")

entrylabel_SAMPLESIZE <- tklabel(f4, text="Sample size")
entrylabel_SIZE <- tklabel(f4,text="Random sample size")


# Position and display sampling parameters on the grid:
tkgrid(entrylabel_NOSAMP,entrylabel_SAMP, entrylabel_RAND)
tkgrid(entry_NOSAMP, entry_SAMP, entry_RAND)
tkgrid(tklabel(f4,text="    "),entrylabel_SAMPLESIZE, entrylabel_SIZE)
tkgrid(tklabel(f4,text="    "),entry_SAMPLESIZE, entry_SIZE)
tkgrid(tklabel(f4,text="    ")) # blank line for aesthetic purposes

# Tooltips for the above
tk2tip(entrylabel_SAMP, "Specify whether the texts in the corpus should be divided in equal-sized samples.")
tk2tip(entrylabel_SAMPLESIZE, "Specify the size for the samples (expressed in words). \nOnly relevant when normal sampling is switched on.")
tk2tip(entrylabel_RAND, "When the analyzed texts are significantly unequal in length, \nit is not a bad idea to prepare samples as randomly chosen *bags of words*. \nIf this option is switched on, the desired size of a sample should be indicated.")
tk2tip(entrylabel_SIZE, "Specify the random sample size. \nOnly relevant when random sampling is switched on.")
tk2tip(entrylabel_NOSAMP, "No internal sampling will be performed: entire texts are considered as samples.")

# next row: OUTPUT
#
cb_SCRN <- tkcheckbutton(f5)
cb_PDF <- tkcheckbutton(f5)
cb_JPG <- tkcheckbutton(f5)
cb_EMF <- tkcheckbutton(f5)
cb_PNG <- tkcheckbutton(f5)
cb_PLOT.RESET <- tkcheckbutton(f5)
entry_COLOR <- tkradiobutton(f5)
entry_GRAY <- tkradiobutton(f5)
entry_BW <- tkradiobutton(f5)
entry_LABELS <-tkradiobutton(f5)
entry_POINTS <-tkradiobutton(f5)
entry_BOTH <-tkradiobutton(f5)
entry_CLASSIC <- tkradiobutton(f5)
entry_LOADINGS <- tkradiobutton(f5)
entry_TECHNICAL <- tkradiobutton(f5)
entry_SYMBOLS <- tkradiobutton(f5)
cb_TITLE <- tkcheckbutton(f5)
cb_HORIZ <- tkcheckbutton(f5)
cb_TABLESAVE <- tkcheckbutton(f5)
cb_FEATURESAVE <- tkcheckbutton(f5)
cb_FREQSAVE <-tkcheckbutton(f5)

tkconfigure(entry_LABELS, variable=text.id.on.graphs, value="labels")
tkconfigure(entry_POINTS, variable=text.id.on.graphs, value="points")
tkconfigure(entry_BOTH, variable=text.id.on.graphs, value="both")
entry_MARGINS <- tkentry(f5,textvariable=add.to.margins,width="8")
entry_OFFSET <- tkentry(f5,textvariable=label.offset,width="8")
entry_PLOT.HEIGHT <- tkentry(f5,textvariable=plot.custom.height,width="8")
entry_PLOT.WIDTH <- tkentry(f5,textvariable=plot.custom.width,width="8")
entry_PLOT.FONT <- tkentry(f5,textvariable=plot.font.size,width="8")
entry_PLOT.LINE <- tkentry(f5,textvariable=plot.line.thickness,width="8")
tkconfigure(cb_SCRN,variable=display.on.screen)
tkconfigure(cb_PDF,variable=write.pdf.file)
tkconfigure(cb_JPG,variable=write.jpg.file)
tkconfigure(cb_EMF,variable=write.emf.file)
tkconfigure(cb_PNG,variable=write.png.file)
tkconfigure(entry_COLOR, variable=colors.on.graphs, value="colors")
tkconfigure(entry_GRAY, variable=colors.on.graphs, value="grayscale")
tkconfigure(entry_BW, variable=colors.on.graphs, value="black")
tkconfigure(cb_PLOT.RESET,variable=plot.options.reset)
tkconfigure(entry_CLASSIC, variable=pca.visual.flavour, value="classic")
tkconfigure(entry_LOADINGS, variable=pca.visual.flavour, value="loadings")
tkconfigure(entry_TECHNICAL, variable=pca.visual.flavour, value="technical")
tkconfigure(entry_SYMBOLS, variable=pca.visual.flavour, value="symbols")
tkconfigure(cb_TITLE,variable=titles.on.graphs)
tkconfigure(cb_HORIZ,variable=dendrogram.layout.horizontal)
tkconfigure(cb_TABLESAVE,variable=save.distance.tables)
tkconfigure(cb_FEATURESAVE,variable=save.analyzed.features)
tkconfigure(cb_FREQSAVE,variable=save.analyzed.freqs)


#
cblabel_SCRN <- tklabel(f5, text="     Onscreen     ")
cblabel_PDF <- tklabel(f5,text="       PDF        ")
cblabel_JPG <- tklabel(f5,text="       JPG        ")
cblabel_EMF <- tklabel(f5,text="       EMF        ")
cblabel_PNG <- tklabel(f5,text="       PNG        ")
entrylabel_COLOR <- tklabel(f5,text="Colors")
entrylabel_GRAY <- tklabel(f5,text="Grayscale")
entrylabel_BW <- tklabel(f5,text="Black")
cblabel_PLOT.RESET <- tklabel(f5,text="Set default")
entrylabel_PLOT.HEIGHT <- tklabel(f5,text="Plot height")
entrylabel_PLOT.WIDTH <- tklabel(f5,text="Plot width")
entrylabel_PLOT.FONT <- tklabel(f5,text="Font size")
entrylabel_PLOT.LINE <- tklabel(f5,text="Line width")
entrylabel_LABELS <- tklabel(f5,text="Labels")
entrylabel_POINTS <- tklabel(f5,text="Points")
entrylabel_BOTH <- tklabel(f5,text="Both")
entrylabel_MARGINS <- tklabel(f5,text="Margins")
entrylabel_OFFSET <- tklabel(f5,text="Label offset")
entrylabel_CLASSIC <- tklabel(f5,text="      Classic    ")
entrylabel_LOADINGS <- tklabel(f5,text="      Loadings   ")
entrylabel_TECHNICAL <- tklabel(f5,text="      Technical  ")
entrylabel_SYMBOLS <- tklabel(f5,text="      Symbols    ")
cblabel_TITLE <- tklabel(f5,text="      Titles      ")
cblabel_HORIZ <- tklabel(f5,text="Horizontal CA tree")
cblabel_TABLESAVE <- tklabel(f5,text="Save distance table")
cblabel_FEATURESAVE <- tklabel(f5,text="Save features")
cblabel_FREQSAVE <- tklabel(f5,text="Save frequencies")

#
tkgrid(tklabel(f5,text="     GRAPHS:"), cblabel_SCRN,cblabel_PDF, cblabel_JPG,cblabel_EMF,cblabel_PNG,columnspan=5)
tkgrid(tklabel(f5,text="            "), cb_SCRN,cb_PDF,cb_JPG,cb_EMF,cb_PNG,columnspan=5)
tkgrid(tklabel(f5,text="    ")) # blank line for aesthetic purposes
tkgrid(tklabel(f5,text="  PLOT AREA:"), cblabel_PLOT.RESET,entrylabel_PLOT.HEIGHT, entrylabel_PLOT.WIDTH,entrylabel_PLOT.FONT,entrylabel_PLOT.LINE,columnspan=5)
tkgrid(tklabel(f5,text="            "), cb_PLOT.RESET,entry_PLOT.HEIGHT,entry_PLOT.WIDTH,entry_PLOT.FONT,entry_PLOT.LINE,columnspan=5)
tkgrid(tklabel(f5,text="            "), tklabel(f5,text="    "),entrylabel_COLOR,entrylabel_GRAY,entrylabel_BW,cblabel_TITLE,columnspan=5)
tkgrid(tklabel(f5,text="            "), tklabel(f5,text="    "),entry_COLOR,entry_GRAY,entry_BW,cb_TITLE,columnspan=5)
tkgrid(tklabel(f5,text="    ")) # blank line for aesthetic purposes
tkgrid(tklabel(f5,text="    PCA/MDS:"), entrylabel_LABELS, entrylabel_POINTS, entrylabel_BOTH, entrylabel_MARGINS,entrylabel_OFFSET, columnspan=5)
tkgrid(tklabel(f5,text="            "), entry_LABELS, entry_POINTS, entry_BOTH, entry_MARGINS,entry_OFFSET, columnspan=5)
tkgrid(tklabel(f5,text="    ")) # blank line for aesthetic purposes
tkgrid(tklabel(f5,text="PCA FLAVOUR:"), entrylabel_CLASSIC, entrylabel_LOADINGS, entrylabel_TECHNICAL, entrylabel_SYMBOLS, columnspan=5)
tkgrid(tklabel(f5,text="            "), entry_CLASSIC, entry_LOADINGS, entry_TECHNICAL, entry_SYMBOLS, columnspan=5)
tkgrid(tklabel(f5,text="    VARIOUS:"), cblabel_HORIZ,cblabel_TABLESAVE,cblabel_FEATURESAVE,cblabel_FREQSAVE,columnspan=5)
tkgrid(tklabel(f5,text="            "), cb_HORIZ,cb_TABLESAVE,cb_FEATURESAVE,cb_FREQSAVE,columnspan=5)
tkgrid(tklabel(f5,text="    ")) # blank line for aesthetic purposes

# Tooltips for the above
tk2tip(cblabel_SCRN, "Select to have your diagram(s) displayed on R's standard graphics device.")
tk2tip(cblabel_PDF, "Select to save your diagram(s) as (a) PDF file(s).")
tk2tip(cblabel_JPG, "Select to save your diagram(s) as (a) JPG file(s).")
tk2tip(cblabel_EMF, "Select to save your diagram(s) as (a) EMF file(s). \nOnly works in Windows.")
tk2tip(cblabel_PNG, "Select to save your diagram(s) as (a) PNG file(s). \nProbably the best option for quality.")
tk2tip(cblabel_PLOT.RESET, "Restore graphic parameters to safe default values \n(7x7 inches, 10 points font size, normal line width).")
tk2tip(entrylabel_PLOT.HEIGHT, "Set custom plot height (in inches).")
tk2tip(entrylabel_PLOT.WIDTH, "Set custom plot width (in inches).")
tk2tip(entrylabel_PLOT.FONT, "Set custom font size (in points).")
tk2tip(entrylabel_PLOT.LINE, "Set custom line width \n(in R units, 1 is default).")
tk2tip(entrylabel_COLOR, "Select to have automatic color coding for tip labels by author.")
tk2tip(entrylabel_GRAY, "Select to have automatic color coding (in grayscale) for tip labels by author.")
tk2tip(entrylabel_BW, "Select to have a black & white graph.")
tk2tip(cblabel_TITLE, "Select to have automatic titles (folder name, analysis type, \nMFW settings, analysis type etc.) on your graph(s).")
tk2tip(entrylabel_LABELS, "Use only labels on MDS/PCA plots.")
tk2tip(entrylabel_POINTS, "Use only points on MDS/PCA plots.")
tk2tip(entrylabel_BOTH, "Use both labels and points on MDS/PCA plots.")
tk2tip(entrylabel_MARGINS, "Set custom margin size \n(in percentage of plot area).")
tk2tip(entrylabel_OFFSET, "Set custom offset between label and point \n(in percentage of plot area.")
tk2tip(entrylabel_CLASSIC, "Original PCA visualization using (colored) sample names") 
tk2tip(entrylabel_LOADINGS, "Display PCA feature (word etc.) loadings.") 
tk2tip(entrylabel_TECHNICAL, "Technical grayscale PCA visualization, showing feature loadings as well as a PC barplot.\nPotentially useful for greyscale printing in traditional publications.")
tk2tip(entrylabel_SYMBOLS, "Select to display the samples in your PCA with a group symbol (instead of their entire name).\n Potentially useful when dealing with lots of samples.")
tk2tip(cblabel_HORIZ, "Select to have your Cluster Analysis graph oriented horizontally. \nProbably the better option for clarity.")
tk2tip(cblabel_TABLESAVE, "Save final distance table(s) in separate text file(s).")
tk2tip(cblabel_FEATURESAVE, "Save final feature (word, n-gram) list(s), e.g. the words actually used in the analysis.")
tk2tip(cblabel_FREQSAVE, "Save frequency table(s) in separate text file(s).")

# next row: the OK button
#
# button_1 <- tkbutton(tt,text="     OK     ",command=push_OK,relief="groove")
# tkbind(button_1,"<Return>",push_OK) 
# tkgrid(button_1,columnspan="10")
tkgrid(tklabel(tt,text="    ")) # blank line (i.e., bottom margin)


##########

repeat{
  if(cancel_pause){
    analyzed.features <- as.character(tclvalue(analyzed.features))
    ngram.size <- as.numeric(tclvalue(ngram.size))
    corpus.format <- as.character(tclvalue(corpus.format))
    mfw.min <- as.numeric(tclvalue(mfw.min))
    mfw.max <- as.numeric(tclvalue(mfw.max))
    mfw.incr <- as.numeric(tclvalue(mfw.incr))
    start.at <- as.numeric(tclvalue(start.at))
    culling.min <- as.numeric(tclvalue(culling.min))
    culling.max <- as.numeric(tclvalue(culling.max))
    culling.incr <- as.numeric(tclvalue(culling.incr))
    use.existing.freq.tables <- as.logical(as.numeric(tclvalue(use.existing.freq.tables)))
    use.existing.wordlist <- as.logical(as.numeric(tclvalue(use.existing.wordlist)))
    interactive.files <- as.logical(as.numeric(tclvalue(interactive.files)))
    use.custom.list.of.files <- as.logical(as.numeric(tclvalue(use.custom.list.of.files)))
    analysis.type <- as.character(tclvalue(analysis.type))
    delete.pronouns <- as.logical(as.numeric(tclvalue(delete.pronouns)))
    display.on.screen <- as.logical(as.numeric(tclvalue(display.on.screen)))
    write.pdf.file <- as.logical(as.numeric(tclvalue(write.pdf.file)))
    write.jpg.file <- as.logical(as.numeric(tclvalue(write.jpg.file)))
    write.emf.file <- as.logical(as.numeric(tclvalue(write.emf.file)))
    write.png.file <- as.logical(as.numeric(tclvalue(write.png.file)))
    colors.on.graphs <- as.character(tclvalue(colors.on.graphs))
    pca.visual.flavour <- as.character(tclvalue(pca.visual.flavour))
    titles.on.graphs <- as.logical(as.numeric(tclvalue(titles.on.graphs)))
    dendrogram.layout.horizontal <- as.logical(as.numeric(tclvalue(dendrogram.layout.horizontal)))
    save.distance.tables <- as.logical(as.numeric(tclvalue(save.distance.tables)))
    save.analyzed.features <- as.logical(as.numeric(tclvalue(save.analyzed.features)))
    save.analyzed.freqs <- as.logical(as.numeric(tclvalue(save.analyzed.freqs)))
    sampling <- as.character(tclvalue(sampling))
    sample.size <- as.numeric(tclvalue(sample.size))
    length.of.random.sample <- as.numeric(tclvalue(length.of.random.sample))
    mfw.list.cutoff <- as.numeric(tclvalue(mfw.list.cutoff))
    distance.measure <- as.character(tclvalue(distance.measure))
    corpus.lang <- as.character(tclvalue(corpus.lang))
    consensus.strength <- as.numeric(tclvalue(consensus.strength))
    plot.options.reset <- as.logical(as.numeric(tclvalue(plot.options.reset)))
    plot.custom.height <- as.numeric(tclvalue(plot.custom.height))
    plot.custom.width <- as.numeric(tclvalue(plot.custom.width))
    plot.font.size <- as.numeric(tclvalue(plot.font.size))
    plot.line.thickness <- as.numeric(tclvalue(plot.line.thickness))
    text.id.on.graphs <- as.character(tclvalue(text.id.on.graphs))
    add.to.margins <- as.numeric(tclvalue(add.to.margins))
    label.offset <- as.numeric(tclvalue(label.offset))
  break
  }
.Tcl("font delete myDefaultFont")
}

} # <-- here the option "interactive.mode.with.GUI == TRUE" is completed

# #################################################
# GUI module explicit feliciter (Phew!)
# #################################################



# #############################################################################
# Final settings (you are advised rather not to change them)
# #############################################################################


# The chosen language option should be followed by an assignment of 
# the appropriate set of pronouns. The following code is responsible for it

  if(corpus.lang == "English")
      pronouns = eng.pronouns
  if(corpus.lang == "Polish")
      pronouns = pol.pronouns 
  if(corpus.lang == "Latin")
      pronouns = lat.pronouns
  if(corpus.lang == "Latin.corr")
    pronouns = lat.pronouns
  if(corpus.lang == "French")
      pronouns = fra.pronouns
  if(corpus.lang == "German" )
      pronouns = ger.pronouns
  if(corpus.lang == "Italian")
      pronouns = ita.pronouns
  if(corpus.lang == "Hungarian")
      pronouns = hun.pronouns
  if(corpus.lang == "Dutch")
    pronouns = dut.pronouns

# Since it it not so easy to perform, say, 17.9 iterations, or analyze
# 543.3 words, the code below rounds off all numerical variables to 
# the nearest positive integers, to prevent you from making silly jokes 
# with funny settings. (OK, it is still possible to crash the script in 
# more ways than one, but you will have to find them on your own).

  mfw.min = round(mfw.min)
  mfw.max = round(mfw.max)
  mfw.incr = round(mfw.incr)
  start.at = round(start.at)
  culling.min = round(culling.min)
  culling.max = round(culling.max)
  culling.incr = round(culling.incr)
  mfw.list.cutoff = round(mfw.list.cutoff)
  sample.size = round(sample.size)


# resetting the default plot area (if an appropriate option has been chosen)
if(plot.options.reset == TRUE) {
  plot.custom.height = 7
  plot.custom.width = 7
  plot.font.size = 10
  plot.line.thickness = 1
  plot.options.reset = FALSE
  }



# If TXM compatibility mode has been chosen, other options need to be switched off
if(txm.compatibility.mode == TRUE) {
  # checking if a frequency table has been passed from TXM to R
  if(exists("txm.generated.freq.table") == TRUE) {
    # inheriting the table from TXM
    frequencies.0.culling = t(variable.name)
    # transposing the table
    frequencies.0.culling = frequencies.0.culling[-1,]
    # set the variable use.existing.freq.tables to skip uploading corpus files
    use.existing.freq.tables == TRUE
  } else {
     cat("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
    "Oops! To use TXM compatibility mode, you have to launch TXM first!\n",
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
    stop("Incorrect input data")
  }
}



# Finally, we want to save some of the variable values for later use;
# they are automatically loaded into the GUI at the next run of the script.
cat("",file="config.txt",append=F)
var.name<-function(x) { 
      if(is.character(x)==TRUE) {
      cat(paste(deparse(substitute(x))," = \"",x,"\"", sep=""),file="config.txt",sep="\n",append=T)
        } else {
          cat(paste(deparse(substitute(x)),x, sep=" = "),file="config.txt",sep="\n",append=T) }
        } 
var.name(corpus.format)
var.name(corpus.lang)
var.name(analyzed.features)
var.name(ngram.size)
var.name(mfw.min)
var.name(mfw.max)
var.name(mfw.incr)
var.name(start.at)
var.name(culling.min)
var.name(culling.max)
var.name(culling.incr)
var.name(mfw.list.cutoff)
var.name(delete.pronouns)
var.name(analysis.type)
var.name(use.existing.freq.tables)
var.name(use.existing.wordlist)
var.name(use.custom.list.of.files)
var.name(consensus.strength)
var.name(distance.measure)
var.name(display.on.screen)
var.name(write.pdf.file)
var.name(write.jpg.file)
var.name(write.emf.file)
var.name(write.png.file)
var.name(save.distance.tables)
var.name(save.analyzed.features)
var.name(save.analyzed.freqs)
var.name(colors.on.graphs)
var.name(titles.on.graphs)
var.name(dendrogram.layout.horizontal)
var.name(pca.visual.flavour)
var.name(sampling)
var.name(sample.size)
var.name(length.of.random.sample)
var.name(sampling.with.replacement)
var.name(plot.custom.height)
var.name(plot.custom.width)
var.name(plot.font.size)
var.name(plot.line.thickness)
var.name(label.offset)
var.name(add.to.margins)
var.name(text.id.on.graphs)


# #############################################################################


# no additional margin will be added on PCA/MDS plots, unless points & labels
# are switched to be shown together
if(text.id.on.graphs != "both") {
  label.offset = 0
  }


# if a chosen plot area is really large (and a bitmap output has been chosen),
# a warning will appear
if(write.jpg.file == TRUE || write.emf.file == TRUE || write.png.file == TRUE){
  # if the desired height*width (at 300 dpi) exceeds 36Mpx
  if(300*plot.custom.width * 300*plot.custom.height > 36000000) {
    cat("\nYou have chosen a bitmap output format and quite a large plot area\n")
    cat("of", plot.custom.width, "by", plot.custom.height, "inches. Producing some",
        as.integer(300*plot.custom.width * 300*plot.custom.height / 1000000),
        "Megapixels will take a good while.\n\n")
    cat("  i - ignore this warning and continue with the current settings\n")
    cat("  p - use pdf format instead of a bitmap (default)\n")
    cat("  s - shrink the plot area to a reasonable size of 20x20 inches\n")
    cat("  a - abort the script\n")
    # reading from the prompt
    answer = readline("\n[i/p/s/a]  ")
    if(tolower(answer) == "a") {
      stop("The script stopped by the user")
    } else if(tolower(answer) == "i") {
      cat("Okay (but honestly, do you really need such a large plot?)\n")
    } else if(tolower(answer) == "s") {
      cat("The plot area will be shrunken to 20x20 inches\n")
      plot.custom.width = 20
      plot.custom.height = 20
    } else {
      cat("Withdrawing from the bitmap output, performing pdf instead\n")
      write.jpg.file = FALSE
      write.emf.file = FALSE
      write.png.file = FALSE
      write.pdf.file = TRUE
    }
  }
}

# #############################################################################



 





# #################################################
# FUNCTIONS:


# #################################################
# Function for combining single features (words
# or characters) into n-grams, or strings of n elements;
# e.g. character 2-grams of the sentence "This is a sentence"
# are as follows: "th", "hi", "is", "s ", " i", "is", etc.
# Required argument: name of the vector of words/chars
# #################################################
make.ngrams = function(input.text) {
  txt = c()
  if(ngram.size > 1) {
    txt = input.text
    for(n in 2:ngram.size) {
    txt = paste(txt[1:(length(txt)-1)],input.text[n:length(input.text)])
    }
  } else {
  # if n-gram size is set to 1, then nothing will happen
  txt = input.text
  }
return(txt)
}


# #################################################
# The generic function for splitting a given input text into
# single words (chains of characters delimited with
# spaces or punctuation marks). Alternatively, 
# you can replace it with another rule.
# Required argument: name of the text to be split.
# ATTENTION: this is (almost) the only piece of coding in this script
# that dependens on the operating system used
# #################################################
split.into.words = function(input.text) {
  # splitting into units specified by regular expression; here, 
  # all sequences between non-letter characters are assumed to be words:
  if(Sys.info()[["sysname"]] == "Windows") { 
    ### Windows
    tokenized.text = c(unlist(strsplit(input.text, "\\W+|_+",perl=T)))
  } else {
    ### Linux, Mac
    tokenized.text = c(unlist(strsplit(input.text, "[^[:alpha:]]+")))
  }
  return(tokenized.text)
}


# #################################################
# Function that carries out the necessary modifications
# for feature selection: convert a sample into
# the type of sequence needed (ngrams etc.) and
# returns the new list of items
# Argument: a vector of words (or chars)
# #################################################
sample.to.features = function(sample){
  # 1. for splitting a given input text into
  # single words (chains of characters delimited with
  # spaces or punctuation marks). Alternatively, 
  # you can replace it with another rule.
  #
  # Splitting the sample into chars (if analyzed.features was set to "c")
  if(analyzed.features == "c") {
    sample = paste(sample, collapse=" ")
    sample = unlist(strsplit(sample,""))
  }
  # 2. making n-grams (if an appropriate option has been chosen):
  if(ngram.size > 1) {
    sample = make.ngrams(sample)
  }
  return(sample)
}


# #################################################
# Function for splitting a given input text into
# single words (chains of characters delimited with
# spaces or punctuation marks). There is also an option
# of splitting the text into characters and/or performing
# splitting into n-grams (see above)
# #################################################
parse.text = function(input.text) {
  # loading the file; optionally, fiddling with apostrophes and contractions:
  #
  # this is the standard procedure of splitting input texts
  if(corpus.lang != "English.contr" && corpus.lang != "English.all") {
    tokenized.text = split.into.words(input.text)
  }
  # if the Latin option with adjusting the v/u letters is on,
  # this smashes the distinction and converts both types to the letter u
  if(corpus.lang == "Latin.corr") {
    tokenized.text = gsub("v","u",tokenized.text)
  }
  # this code is used for English corpora only
  if(corpus.lang == "English.contr" || corpus.lang == "English.all") {
    # replacing non-ASCII apostrophes with simple ' (standard ASCII char)
    tokenized.text = gsub("’","'",input.text)
    # getting rid of contractions ('t, 's, 've, 'd, 'll, 'em, 'im) by replacing
    # their apostrophes with ^ (other apostrophes will not be replaced);
    # Of course, if your corpus is Cockney, you should edit the 
    # "([tsdm]|ll|ve|em|im)" statement accordingly.
    tokenized.text = gsub("([[:alpha:]])'([tsdm]|ll|ve|em|im)\\b","\\1^\\2",
                            tokenized.text) #'
    # adding spaces around dashes (to distinguish dashes and hyphens)
    tokenized.text = gsub("[-]{2,5}"," -- ",tokenized.text)
    # depending on which option was swithed on, either the contractions are
    # kept, or all the peculiarities, i.e. both contractions and hyphens
    if(corpus.lang == "English.contr") {
      tokenized.text=c(unlist(strsplit(tokenized.text,"[^[:alpha:]^]+")))
  }
    if(corpus.lang == "English.all") {
      tokenized.text=c(unlist(strsplit(tokenized.text,"[^[:alpha:]^-]+")))
      # trying to clean the remaining dashes:
      tokenized.text = gsub("^[-]+$","",tokenized.text)
  }
  }
  # trying to avoid empty strings:
  tokenized.text = tokenized.text[nchar(tokenized.text)>0]
  # trying to get rid of non-letter characters:
  tokenized.text = tokenized.text[grep("[^[:digit:]]",tokenized.text)]
  # sanity check for text length: abort if the current text is extremely
  # short or at least shorter than the specified sample size
  if (length(tokenized.text) < 10 || 
      (sampling == "normal.sampling" && length(tokenized.text) < sample.size) || 
      (sampling == "random.sampling" && length(tokenized.text) < sample.size)) {
    cat("\n\n",file, "\t", "This sample is too short!", "\n\n")
    setwd(".."); stop("Corpus error...")
  }
  # at this point, each text in the corpus has been tokenized
  # into an array of tokens which we can divide into samples
  samples.from.text = list()
  if (sampling == "normal.sampling"){
    # initialize variables to sample the text
    text.length = length(tokenized.text)
    number.of.samples = floor(text.length/(sample.size))
    cat(paste("\t", "- text length (in words): ", text.length, "\n", sep=""))
    cat(paste("\t", "- nr. of samples: ", number.of.samples, "\n", sep=""))
    cat(paste("\t", "- nr. of words dropped at the end of the file: ", text.length-(number.of.samples*sample.size), "\n", sep=""))
    # iterate over the samples:
    current.start.index = 1
    for(sample.index in 1:number.of.samples) {
      current.sample = tokenized.text[current.start.index:(current.start.index+sample.size-1)]
      current.sample = sample.to.features(current.sample)
      # flush current sample:
      samples.from.text[[sample.index]] = current.sample
      # increment index for next iteration
      current.start.index = current.start.index + sample.size
      current.sample = c()
    }
  } else if(sampling == "random.sampling"){
    # if random sampling was chosen, the text will be randomized and a sample of a given length will be excerpted
    current.sample = head(sample(tokenized.text, replace = sampling.with.replacement), length.of.random.sample)
    current.sample = sample.to.features(current.sample)
    samples.from.text[[1]] = current.sample 
  } else if (sampling == "no.sampling"){
    # entire texts will be used as a sample (regardless of its length)
    current.sample = tokenized.text
    current.sample = sample.to.features(current.sample)
    samples.from.text[[1]] = current.sample
 }
  return(samples.from.text)
}


# #################################################
# Function for adjusting different input formats:
# xml (TEI) in two variants, html, and plain text files.
# Required argument: name of the text to pre-process
# #################################################
delete.markup = function(input.text) {
  if(corpus.format == "xml" || corpus.format == "xml.drama") {
    # getting rid of the TEI header (if it exists)
    if(length(grep("</teiheader>",input.text)) > 0) {
      input.text = input.text[-c(1:(grep("</teiheader>",input.text)))]
      }
    # the whole text into one (very) long line
    preprocessed.text = paste(input.text, collapse=" ")
    # getting rid of dramatis personae
    if(corpus.format == "xml.drama"){
      preprocessed.text = gsub("<speaker>.*?</speaker>","",preprocessed.text)
      }
    # getting rid of comments and (editorial) notes
    preprocessed.text = gsub("<note.*?</note>","",preprocessed.text)
    # getting rid of all the remaining tags
    preprocessed.text = gsub("<.*?>","",preprocessed.text)
  }
  if(corpus.format == "html") {
    # getting rid of html header (if exists)
    if(length(grep("<body",input.text)) > 0) {
      input.text = input.text[-c(1:(grep("<body",input.text)))]
      }
    # the whole text into one (very) long line
    preprocessed.text = paste(input.text, collapse=" ")
    # getting rid of links (menus and similar stuff should be deleted, hopefully)
    preprocessed.text = gsub("<a href.*?/a>","",preprocessed.text)
    # getting rid of all the remaining tags
    preprocessed.text = gsub("<.*?>","",preprocessed.text)
  } else {
  preprocessed.text = input.text
  }
return(preprocessed.text)
}



# #################################################
# Function for graph auto-coloring; depending on
# the user's choice, it assigns colors or grayscale tones
# to matching strings of characters in texts' names
# (as a delimiter, the underscore character is used);
# alternatively, all the labels can be marked black. 
# Required argument: a vector of text labels
# Optional argument: col="colors" || "grayscale" || "back"
# #################################################
assign.plot.colors = function(names.of.the.texts,col="colors") {
  if(col == "black") {
    colors.of.pca.graph = "black"
  } else {
    color.numeric.values = c(1)
    current.color = 1
    # a loop for matching the subsequent strings of chars
    for(w in 2:length(names.of.the.texts)) {
      # if two samples have the same id, don't change the color
      if(gsub("_.*","",names.of.the.texts)[w] %in%  
                       gsub("_.*","",names.of.the.texts[1:(w-1)]) == TRUE) {
        find.color = which(gsub("_.*","",names.of.the.texts) == 
                               gsub("_.*","",names.of.the.texts)[w])[1]
        current.color = color.numeric.values[find.color]
      # if the samples differ, assign the next color (actually, the color's id)
      } else {
        current.color = max(color.numeric.values) + 1
      }
    # append the recent color to the final vector of color values
    color.numeric.values = c(color.numeric.values, current.color)
    }
  # define a vector of available colors, if an appropriate option was chosen
  if(col == "colors") {
    available.colors = rep(c("red","green","blue","black","orange","purple",
      "darkgrey","brown","maroon4","mediumturquoise","gold4", "deepskyblue",
      "yellowgreen","grey","chartreuse4", "khaki", "navy", "palevioletred",
      "greenyellow", "darkolivegreen4", "chocolate4"),10)
    }
  # define a vector of gray tones, instead of colors
  if(col == "grayscale") {
    number.of.colors.required = max(color.numeric.values)
    available.colors = gray(seq(0,0.7,0.7/(number.of.colors.required-1)))
  }
  # produce the final vector of colors (or gray tones)
  colors.of.pca.graph = available.colors[c(color.numeric.values)]
  }
return(colors.of.pca.graph)
}




# #################################################
# Function that finds out the coordinates
# of scatterplots; it computes the extreme x and y
# values, adds some margins, and optionally extends
# the top margin if a plot uses sample labels
# Required arguments: (1) a vector of x coordinates,
# (2) a vector of y coordinates, optionally with names;
# optional arguments: (3) additional margins (expressed 
# in % of actual plot area), (4) label offset (in %)
# #################################################
define.plot.area = function(x.coord,y.coord,xymargins=2,v.offset=0) {
  # get horizontal plotting area (reasonable margins added on both sides):
  # (1) how long are the extreme samples' names
  left.label.length = nchar(names(x.coord)[(order(x.coord)[1])])
  right.label.length = nchar(names(x.coord)[(order(x.coord,decreasing=T)[1])])
  # (2) checking if the sample labels really exist
    if(length(left.label.length) == 0) {
      left.label.length = 0}
    if(length(right.label.length) == 0) {
      right.label.length = 0}
  # (3) x axis expansion (0.5% for each character of the extreme samples' names)
  left.expansion = left.label.length * 0.005
  right.expansion = right.label.length * 0.005
  # (4) size of the x axis
  x.axis.size = abs(max(x.coord) - min(x.coord))
  # (5) finally, get both x coordinates
  min.x = min(x.coord) - (left.expansion + 0.01 * xymargins) * x.axis.size
  max.x = max(x.coord) + (right.expansion + 0.01 * xymargins) * x.axis.size 
  #
  # get vertical plotting area (plus optional top margin):
  # (1) size of the y axis
  y.axis.size = abs(max(y.coord) - min(y.coord))
  # (2) top margin (added to fit the samples' labels)
  top.offset = 0.005 * y.axis.size
  # (3) finally, get both y coordinates
  min.y = min(y.coord) - 0.01 * xymargins * y.axis.size
  max.y = max(y.coord) + (0.01 * label.offset + 0.01 * xymargins) * y.axis.size
  #
  plot.area = list(c(min.x, max.x), c(min.y, max.y))
return(plot.area)
}


# #################################################








# #################################################
# the module for loading a corpus from the text files;
# it can be omitted if the frequency table already exists
# (then "use.existing.freq.tables" should be set 
# to TRUE in the preamble of the script/GUI)
# #################################################
#
# Checking: (1) whether to produce a new frequency table or to use 
# the existing one; (2) whether the tables are stored in memory or 
# written into files.
# If you have chosen using the existing table and it does not exist,
# available, then your choice will be ignored and the table will be 
# created from scratch.


# checking if an appropriate frequency table exists
if(exists("frequencies.0.culling") == FALSE 
              && file.exists("table_with_frequencies.txt") == FALSE ) {
  use.existing.freq.tables = FALSE
}


if(use.existing.freq.tables == TRUE) { 
      if(exists("frequencies.0.culling")) {
      cat("\n", "using frequency table stored as variables...", "\n")
        } else {
          cat("\n", "reading file with frequency table...", "\n")
          frequencies.0.culling = t(read.table("table_with_frequencies.txt"))
          cat("\n", "frequency table loaded successfully", "\n\n")
        }
      # extracting names of the texts
      corpus.filenames = rownames(frequencies.0.culling)
      #
      # checking whether an existing wordlist should be used
      if (use.existing.wordlist == TRUE && file.exists("wordlist.txt") == TRUE){
          cat("\n", "reading a wordlist from file...", "\n")
          mfw.list.of.all = scan("wordlist.txt",what="char",sep="\n")
          mfw.list.of.all = c(grep("^[^#]",mfw.list.of.all,value=TRUE))
          #
          # adjusting the size of the frequency table according to the existing wordlist
          frequencies.0.culling = 
                       frequencies.0.culling[,colnames(frequencies.0.culling) 
                       %in% mfw.list.of.all]
      } else {
          # the wordlist will be created from the existing frequency tables
          mfw.list.of.all = colnames(frequencies.0.culling)
          # some comments into the file containing the wordlist
          cat("# This file contains the words that were used in the table",
          "# of frequencies uploaded from an external file. The current list",
          "# can be used for the next tasks, and for this purpose it can be",
          "# manually revised, edited, deleted, culled, etc.", 
          "# You can either delete unwanted words, or mark them with \"#\"",
          "# -----------------------------------------------------------------",
          "",
          file="wordlist.txt", sep="\n")
          # the current wordlist into a file
          cat(mfw.list.of.all, file="wordlist.txt", sep="\n",append=F)
        }
# if the existing table will not be used, then begin producing the new table
  } else {
#
# Retrieving the names of texts
#
# first, it's possible to choose the files manually
if (interactive.files == TRUE) {
    setwd("corpus")
  corpus.filenames = basename(tk_choose.files(default = "", caption = "Select at least 2 files", multi = TRUE))
  setwd("..")
} else {
  # alternatively, one can use the files listed in "files_to_analyze.txt";
  # the listed files can be separated by spaces, tabs, or newlines
  if(use.custom.list.of.files == TRUE 
      && file.exists("files_to_analyze.txt") == TRUE) { 
    # a message on the screen
    cat("\n")
    cat("external list of files will be used for uploading the corpus\n\n")
    # retrieving the filenames from a file
    corpus.filenames = scan("files_to_analyze.txt",what="char",sep="\n",quiet=T)
    # getting rid of spaces and/or tabs
    corpus.filenames = unlist(strsplit(corpus.filenames,"[ \t]+"))
      # checking whether all the files indicated on the list really exist
      if( length(setdiff(corpus.filenames,list.files("corpus"))) > 0 ){
        # if not, then sent a message and list the suspicious filenames
        cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        cat("the following files have not been found:\n")
        cat(setdiff(corpus.filenames, list.files("corpus")),"\n\n")
        cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
        # use only those files that match
        corpus.filenames = intersect(corpus.filenames, list.files("corpus"))
      }
  } else {
  corpus.filenames = list.files("corpus")
  }
}
#
# Checking whether the required files and subdirectory exist
if(file.exists("corpus")==FALSE) {
    cat("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
    "Hey! The working directory should contain the subdirectory \"corpus\"\n",
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
    stop("Corpus prepared incorrectly")
    }
if(length(corpus.filenames) < 2 && sampling !="normal.sampling")  {
    cat("\n\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n",
    "Ho! The subdirectory \"corpus\" should contain at least two text samples!\n",
    "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n")
    stop("Corpus prepared incorrectly")
    }
#
# loading the corpus from individual text files
loaded.corpus = list()
if (sampling == "normal.sampling"){
  cat(paste("Performing sampling (using sample size = ", sample.size," words)\n", sep=""))
} else if(sampling == "random.sampling"){
  cat(paste("Performing random sampling (using random sample size = ", " words)\n", sep=""))
} else if (sampling == "no.sampling"){
cat(paste("Performing no sampling (using entire text as sample)", "\n", sep=""))
} else {
  stop("Exception raised: something is wrong with the sampling parameter you have specified...")
}
    
# temporarily change the working directory to the corpus directory
setwd("corpus")
for (file in corpus.filenames) {
  cat(paste("Loading ", file, "\t", "...", "\n", sep=""))
  # loading the next file from the list "corpus.filenames"
  current.file = tolower(scan(file,what="char",sep="\n", quiet=T))
  # delete xml/html markup (if applicable)
  current.file = delete.markup(current.file)
  # deleting punctuation, splitting into words:
  samples.from.text = parse.text(current.file)
  # appending the current text to the virtual corpus
  if (sampling == "normal.sampling"){
    for (index in 1:length(samples.from.text)){
      # add the parsed sample to the corpus (and remove the filename extension)
      loaded.corpus[[paste(gsub("(\\.txt$)||(\\.xml$)||(\\.html$)||(\\.htm$)",
                   "",file), "-", index, sep="")]] = samples.from.text[[index]]
      }
  } else {
      loaded.corpus[[file]] = unlist(samples.from.text)
  }
}
# reset the working directory
setwd("..")
cat(paste("Total nr. of samples in the corpus: ", length(loaded.corpus), "\n"))

#
#
# the directory with corpus must contain enough texts;
# if the number of text samples is lower than 2, the script will abort.
if( (length(corpus.filenames) < 2) & (sampling=="no.sampling") ) {
    cat("\n\n","your corpus folder seems to be empty!", "\n\n")
    stop("corpus error")
}
#
#
# We need a list of the most frequent words used in the current corpus, 
# in descending order, without frequencies (just a list of words). It can be 
# either loaded from a file (then set the option "use.existing.wordlist=TRUE"), 
# or created by the code provided below:
#
    
if (use.existing.wordlist == TRUE && file.exists("wordlist.txt") == TRUE) {
          cat("\n", "reading a wordlist from file...", "\n")
          # loading the wordlist fil  e, changing to lowercase
          mfw.list.of.all = tolower(scan("wordlist.txt",what="char",sep="\n"))
          # getting rid of commented lines in the wordlist file
          mfw.list.of.all = c(grep("^[^#]",mfw.list.of.all,value=TRUE))
} else {
# Extracting all the words used in the corpus
#
wordlist.of.loaded.corpus = c()
  for (file in 1 : length(loaded.corpus)) {
    # loading the next sample from the list "corpus.filenames"
    current.text = loaded.corpus[[file]]
    # putting samples together:
    wordlist.of.loaded.corpus = c(wordlist.of.loaded.corpus, current.text)
#    cat(names(loaded.corpus[file]),"\t","tokenized successfully", "\n")
    }
#
# preparing a sorted frequency list of the whole set
mfw.list.of.all = sort(table(c(wordlist.of.loaded.corpus)),decreasing=T)
  # if the whole list is long, then cut off the tail, as specified in the GUI 
  # by the cutoff value
  if (length(mfw.list.of.all) > mfw.list.cutoff) {
    mfw.list.of.all = mfw.list.of.all[1:mfw.list.cutoff]
  }
# the only thing we need are words ordered by frequency (no frequencies)
mfw.list.of.all = names(mfw.list.of.all)
#
# some comments into the file containing the wordlist
cat("# This file contains the words that were used for building the table",
  "# of frequencies. It can be also used for further tasks, and for this",
  "# purpose it can be manually revised, edited, deleted, culled, etc.", 
  "# You can either delete unwanted words, or mark them with \"#\"",
  "# -----------------------------------------------------------------------",
  "",
      file="wordlist.txt", sep="\n")
# the current wordlist into the "wordlist.txt" file
cat(mfw.list.of.all, file="wordlist.txt", sep="\n",append=T)
#
}   # <----- conditional expr. "use.existing.wordlist" terminates here
#
# blank line on the screen
cat("\n")

#
# #################################################
# FUNCTION: make.parallel.frequency.lists()
# preparing a huge table with all the frequencies (> mwf.list.cutoff).
# Two arguments are required -- a vector with sample names
# and a specified variable where the corpus is stored (in a list)
# #################################################
#
make.parallel.frequency.lists = function(sample.names,current.corpus) {
  freq.list.of.all.the.samples = c()
  freq.list.of.current.sample = c()
    for (sample.name in sample.names) {
    # loading the next sample from the list "sample.names"
    current.sample = current.corpus[[sample.name]]
    # preparing the frequency list of the current text
    raw.freq = table(current.sample) * 100 / length(current.sample)
    # adjusting the frequency list to the main MFW list obtained above
    freq.list.of.current.sample = raw.freq[mfw.list.of.all]
    # taking the names (sc. words) from the main MFW list 
    names(freq.list.of.current.sample) = mfw.list.of.all
    # and inserting the current sample into the general frequency table
    freq.list.of.all.the.samples = rbind(freq.list.of.all.the.samples, freq.list.of.current.sample)
    # a short message on the screen:
#    cat(file, "\t", "excerpted successfully", "\n")
    cat(".")
  }
  # adjusting names of the rows (=samples)
  rownames(freq.list.of.all.the.samples) = c(sample.names)
# the result of the function
return(freq.list.of.all.the.samples)
}
#
# preparing a huge table of all the frequencies for the whole corpus
frequencies.0.culling = make.parallel.frequency.lists(names(loaded.corpus),loaded.corpus)
# all NA values will be adjusted to 0
frequencies.0.culling[which(is.na(frequencies.0.culling))] = 0
#
# getting rid of zero values (this might happen in random sampling
# or when custom wordlist are used)
frequencies.0.culling = frequencies.0.culling[,grep("FALSE",(colSums(frequencies.0.culling))==0)]
#
#
#
# writing the frequency tables to text files (they can be re-used!)
write.table(t(frequencies.0.culling), 
            file="table_with_frequencies.txt", 
            sep="\t",
            row.names=TRUE,
            col.names=TRUE)
}  # <----- conditional expr. "use.existing.freq.tables" terminates here
#
#
# #################################################
# the module for loading the corpus terminates here
# #################################################


# #################################################
# MAIN PROGRAM; the main loop is below
# #################################################

# saving the original mfw.max value in mfw.max.original
# this is useful for graph subtitles
mfw.max.original = mfw.max

# the general counter for various purposes: initialization
number.of.current.iteration = 0

# load the ape library; make an empty bootstrap.results list
# this will be executed only if the bootstrap option is checked
if (analysis.type == "BCT") {
    library(ape)
    bootstrap.list = list()
}


# #################################################
# module for culling (THE MAIN LOOP IN THE PROGRAM)
# #################################################


# testing if desired culling settings are acceptable;
# if too large, it is set to maximum possible
  if(culling.max > 100) {
  culling.max = 100
  }
# if too small, it is set to 0 (i.e. minimal value)
  if(culling.min < 0) {
  culling.min = 0
  }
# avoiding infinite loops
  if(culling.incr <= 1) {
  culling.incr = 10
  }

# #################################################



for(j in (culling.min/culling.incr):(culling.max/culling.incr)) {
current.culling = j * culling.incr

# the beginning of the culling procedure 
raw.list.after.culling = c()

# extracting non-zero values the frequency table.
nonzero.values = frequencies.0.culling > 0


# counting non-zero values
for (y in 1: length(nonzero.values[1,])) {
  raw.list.after.culling = c(raw.list.after.culling, 
              (length(grep("TRUE",nonzero.values[,y])) / 
                     length(nonzero.values[,y])) 
                           >= current.culling/100 
                           )
}
# a raw culling list has no word-identification; let’s change it:
names(raw.list.after.culling) = colnames(frequencies.0.culling)
# a simple sequence of words which have not been culled

list.of.words.after.culling = c(names(raw.list.after.culling[grep("TRUE",raw.list.after.culling)]))

# procedure for deleting pronouns
if (delete.pronouns == TRUE) {
    list.of.words.after.culling = 
      list.of.words.after.culling[!(list.of.words.after.culling %in% pronouns)]
}

# the above list-of-not-culled to be applied to the wordlist:
table.with.all.freqs = frequencies.0.culling[,c(list.of.words.after.culling)]

# the names of the samples are passed to the frequency table
if(use.existing.freq.tables == FALSE) {
  rownames(table.with.all.freqs) = names(loaded.corpus)
}
  
# #################################################
# culling is done, but we are still inside the main loop

# starting the frequency list at frequency rank set in option start.at above
table.with.all.freqs = table.with.all.freqs[,start.at:length(table.with.all.freqs[1,])]


# Testing if the desired MFW number is acceptable,
# if MFW too large, it is set to maximum possible.
  if(mfw.max > length(table.with.all.freqs[1,])) {
  mfw.max = length(table.with.all.freqs[1,])
  }
# if too small, it is set to 1 (i.e., minimal value)
  if(mfw.min < 1) {
  mfw.min = 1
  }
# if culling is too strong, sometimes strange things may happen; let’s block it
  if(mfw.min > mfw.max) {
  mfw.min = mfw.max
  }
# MFW set to mfw.max for a while (it will change later on)
mfw = mfw.max

cat("\n\n")
cat("culling @ ", current.culling,"\t","available words ",mfw.max,"\n")


# #################################################
# z-scores calcutations
# #################################################

if((analysis.type == "CA") || (analysis.type == "BCT") || (analysis.type == "MDS")){
  # calculating z-scores (a message on the screen)
  cat("Calculating z-scores... \n\n")
  # Entropy distance: experimental, but entirely available yet
  # (the results do not really differ than for typical word frequencies)
  #
  #A = t(t(table.with.all.freqs + 1) / colSums(table.with.all.freqs + 1))
  #B =t(t(log(table.with.all.freqs + 2)) / -(colSums(A * log(A))))
  #table.with.all.freqs = B
  #
  # calculating z-scores 
  table.with.all.zscores = scale(table.with.all.freqs)
  table.with.all.zscores = table.with.all.zscores[,]
}

# #################################################
# the internal loop starts here (for i = mfw.min : mfw.max)
# #################################################

# a short message on the screen about distance calculations (when appropriate):
if((analysis.type == "CA") || (analysis.type == "BCT") || (analysis.type == "MDS")){

  if(distance.measure == "CD") {
    cat("Calculating classic Delta distances... \n")
  }
  if(distance.measure == "AL") {
    cat("Calculating Argamon's Delta distances... \n")
  }
  if(distance.measure == "ED") {
    cat("Calculating Eder's Delta distances... \n")
  }
  if(distance.measure == "ES") {
    cat("Calculating Eder's Simple distances... \n")
  }
  if(distance.measure == "MH") {
    cat("Calculating Manhattan distances... \n")
  }
  if(distance.measure == "CB") {
    cat("Calculating Canberra distances... \n")
  }
  if(distance.measure == "EU") {
    cat("Calculating Euclidean distances... \n")
  }
}


cat("MFW used: ")
for(i in (mfw.min/mfw.incr):(mfw.max/mfw.incr)) {
mfw = i * mfw.incr

# for safety reasons, if MFWs > words in samples
if(mfw > length(list.of.words.after.culling) ) {
  mfw = length(list.of.words.after.culling)
}

# the general counter for various purposes
number.of.current.iteration = number.of.current.iteration + 1

# the current task (number of MFW currently analyzed) echoed on the screen
cat(mfw, " ")

# #################################################
# module for calculating distances between texts
# #################################################

if((analysis.type == "CA") || (analysis.type == "BCT") || (analysis.type == "MDS")){
  # calculating Delta distances to a distance matrix
  if(distance.measure == "CD") {
    distance.name.on.graph = "Classic Delta distance"
    distance.name.on.file = "Classic Delta"
    distance.table = 
        as.matrix(dist(table.with.all.zscores[,1:mfw],
        method="manhattan")) / mfw
    }

  # calculating Argamon’s "Linear Delta"
  if(distance.measure == "AL") {
    distance.name.on.graph = "Argamon’s Delta distance"
    distance.name.on.file = "Argamon’s Delta"
    distance.table = 
        as.matrix(dist(table.with.all.zscores[,1:mfw],
        method="euclidean")) / mfw
    }

  # calculating Delta distances with Eder’s modifications
  if(distance.measure == "ED") {
    distance.name.on.graph = "Eder’s Delta distance"
    distance.name.on.file = "Eder’s Delta"
    zscores.plus.e.value = t(t(table.with.all.zscores[,1:mfw])*((1+mfw:1)/mfw))
    distance.table = as.matrix(dist(zscores.plus.e.value,method="manhattan"))
    }

  # calculating Eder’s Simple distance to a distance matrix
  if(distance.measure == "ES") {
    distance.table = 
       as.matrix(dist(sqrt(table.with.all.freqs[,1:mfw]),method="manhattan"))
    distance.name.on.graph = "Eder’s Simple distance"
    distance.name.on.file = "Eder’s Simple"
    }

  # calculating Manhattan distance to a distance matrix
  if(distance.measure == "MH") {
    distance.name.on.graph = "Manhattan distance"
    distance.name.on.file = "Manhattan"
    distance.table = 
         as.matrix(dist(table.with.all.freqs[,1:mfw],method="manhattan"))
    }

  # calculating Canberra distance to a distance matrix
  if(distance.measure == "CB") {
    distance.name.on.graph = "Canberra distance"
    distance.name.on.file = "Canberra"
    distance.table = 
         as.matrix(dist(table.with.all.freqs[,1:mfw],method="canberra"))
    }

  # calculating Euclidean distance to a distance matrix
  if(distance.measure == "EU") {
    distance.name.on.graph = "Euclidean distance"
    distance.name.on.file = "Euclidean"
    distance.table = 
         as.matrix(dist(table.with.all.freqs[,1:mfw],method="euclid"))
    }

  # replaces the names of the samples (the extension ".txt" is cut off)
  rownames(distance.table)=gsub("(\\.txt$)||(\\.xml$)||(\\.html$)||(\\.htm$)",
                        "",rownames(table.with.all.freqs))
  colnames(distance.table)=gsub("(\\.txt$)||(\\.xml$)||(\\.html$)||(\\.htm$)",
                        "",rownames(table.with.all.freqs))
}


# #################################################
# a tiny module for graph auto-coloring: 
# uses the function "assign.plot.colors()"
# #################################################

names.of.the.texts = gsub("(\\.txt)||(\\.xml)||(\\.html)||(\\.htm)","",rownames(table.with.all.freqs))

colors.of.pca.graph = assign.plot.colors(names.of.the.texts,col=colors.on.graphs)





# #################################################
# preparing the graphs
# #################################################

# The name of a given method will appear in the title of the graph
# (if the appropriate option was chosen), and will be pasted into
# a filename of the current job. First, variables are initiated...
name.of.the.method = ""
short.name.of.the.method = ""
mfw.info = mfw 
plot.current.task = function() {NULL}

# getting rid of redundant start.at information
  if(start.at == 1) {
    start.at.info = ""  
    } else {
    start.at.info = paste("Started at",start.at) }
# getting rid of redundant pronoun information
  if(delete.pronouns == TRUE) {
    pronouns.info = paste("Pronouns deleted")  
    } else {
    pronouns.info = "" }
# getting rid of redundant culling information
  if(culling.min == culling.max) {
    culling.info = culling.min 
    } else {
    culling.info = paste(culling.min,"-",culling.max,sep="") }

# prepares a dendrogram for the current MFW value for CA plotting
if(analysis.type == "CA") {
  name.of.the.method = "Cluster Analysis"
  short.name.of.the.method = "CA"
  if(dendrogram.layout.horizontal == TRUE) {
    dendrogram.margins =  c(5,4,4,8)+0.1 
    } else {
    dendrogram.margins = c(8,5,4,4)+0.1 }
  # the following task will be plotted
  plot.current.task = function(){ 
    par(mar=dendrogram.margins)
########################################################################
########################################################################
# color graphs, but using different clustering algorithm (i.e. neighbor joining)
if(nj.cluster.analysis == TRUE) {
  plot(nj(distance.table), font=1, tip.color=colors.of.pca.graph)
  # alternatively, a traditional approach:
} else {
########################################################################
    # clustering the distances stored in the distance.table
    clustered.data = hclust(as.dist(distance.table),"ward")
    # reordering the vector of colors to fit the order of clusters
    colors.on.dendrogram = colors.of.pca.graph[clustered.data$order]
    # converting the clusters into common dendrogram format
    tree.with.clusters = as.dendrogram(clustered.data,hang=0)
    # now, preparing the procedure for changing leaves' color attributes
      # (this snippet is taken from "help(dendrapply)" and slightly adjusted)
      local({
        colLab <<- function(n) {
            if(is.leaf(n)) {
              a <- attributes(n)
              i <<- i+1
              attr(n, "nodePar") <-
                  c(a$nodePar, lab.col = mycols[i], pch = NA)
            }
            n
        }
        mycols = colors.on.dendrogram
        i <- 0
      })
    # adding the attributes to subsequent leaves of the dendrogram,
    # using the above colLab(n) function
    dendrogram.with.colors = dendrapply(tree.with.clusters, colLab)
    # finally, ploting the whole stuff
    plot(dendrogram.with.colors,
           main = graph.title,
           horiz = dendrogram.layout.horizontal) 
    if(dendrogram.layout.horizontal == TRUE) {
      title(sub=graph.subtitle) 
    } else {
      title(sub=graph.subtitle, outer=TRUE, line=-1)  
    }
  }
}}


# prepares a 2-dimensional plot (MDS) for plotting
if(analysis.type == "MDS") {
  name.of.the.method = "Multidimensional Scaling"
  distance.name.on.graph = ""
  distance.name.on.file = ""
  short.name.of.the.method = "MDS" 
  mds.results = cmdscale(distance.table,eig=TRUE)
  # prepare the xy coordinates, add the margins, add the label offset
  xy.coord = mds.results$points[,1:2]
  if(text.id.on.graphs == "both") {
    label.coord = cbind(mds.results$points[,1],(mds.results$points[,2] + (0.01*label.offset*
                      abs(max(mds.results$points[,2]) - min(mds.results$points[,2])))))
    } else {
    label.coord = xy.coord
    }
  plot.area = define.plot.area(mds.results$points[,1],mds.results$points[,2],
                               xymargins=add.to.margins,
                               v.offset=label.offset)
  # define the plotting function needed:
  plot.current.task = function(){ 
    if(text.id.on.graphs == "points" || text.id.on.graphs == "both") {
      plot(xy.coord, type="p", 
           ylab="", xlab="", 
           xlim=plot.area[[1]],ylim=plot.area[[2]],
           main = graph.title,
           sub = graph.subtitle,
           col = colors.of.pca.graph,
           lwd = plot.line.thickness) 
      }
    if(text.id.on.graphs == "labels") {
      plot(xy.coord, type="n", 
           ylab="", xlab="", 
           xlim=plot.area[[1]],ylim=plot.area[[2]],
           main = graph.title,
           sub = graph.subtitle,
           col = colors.of.pca.graph,
           lwd = plot.line.thickness) 
      }
    if(text.id.on.graphs == "labels" || text.id.on.graphs == "both") {
      text(label.coord, rownames(label.coord), col=colors.of.pca.graph) 
      }
    axis(1,lwd=plot.line.thickness)
    axis(2,lwd=plot.line.thickness)
    box(lwd=plot.line.thickness)
  }
}


# prepares Principal Components Analysis (PCA) for plotting
if(analysis.type == "PCV" || analysis.type == "PCR") {
  # set some string information variables
  name.of.the.method = "Principal Components Analysis"
  short.name.of.the.method = "PCA"
  distance.name.on.file = "PCA"
  if(analysis.type == "PCV") {
    pca.results = prcomp(table.with.all.freqs[,1:mfw])
    distance.name.on.graph = "Covariance matrix"
  } else if(analysis.type == "PCR") {
    pca.results = prcomp(table.with.all.freqs[,1:mfw], scale=TRUE)
    distance.name.on.graph = "Correlation matrix"
  }
  # get the variation explained by the PCs:
  expl.var = round(((pca.results$sdev^2)/sum(pca.results$sdev^2)*100),1)
  PC1_lab = paste("PC1 (",expl.var[1],"%)", sep="")
  PC2_lab = paste("PC2 (",expl.var[2],"%)", sep="")

  # prepare the xy coordinates, add the margins, add the label offset
  xy.coord = pca.results$x[,1:2]
  if(text.id.on.graphs == "both") {
    label.coord = cbind(pca.results$x[,1],(pca.results$x[,2] + (0.01*label.offset*
                      abs(max(pca.results$x[,2]) - min(pca.results$x[,2])))))
    } else {
    label.coord = xy.coord
    }
  plot.area = define.plot.area(pca.results$x[,1],pca.results$x[,2],
                               xymargins=add.to.margins,
                               v.offset=label.offset)
  # define the plotting function needed:
  plot.current.task = function(){
    if (pca.visual.flavour == "classic"){
      if(text.id.on.graphs == "points" || text.id.on.graphs == "both") {
        plot(xy.coord,
             type="p",
             xlim=plot.area[[1]],ylim=plot.area[[2]],
             xlab="",ylab=PC2_lab,
             main = graph.title,sub = paste(PC1_lab,"\n",graph.subtitle),
             col=colors.of.pca.graph,
             lwd=plot.line.thickness) 
      }
      if(text.id.on.graphs == "labels") {
        plot(xy.coord,
             type="n",
             xlim=plot.area[[1]],ylim=plot.area[[2]],
             xlab="",ylab=PC2_lab,
             main = graph.title,sub = paste(PC1_lab,"\n",graph.subtitle),
             col=colors.of.pca.graph,
             lwd=plot.line.thickness) 
      }
      abline(h=0, v=0, col = "gray60",lty=2)
      if(text.id.on.graphs == "labels" || text.id.on.graphs == "both") {
        text(label.coord, rownames(pca.results$x), col=colors.of.pca.graph) 
      }
      axis(1,lwd=plot.line.thickness)
      axis(2,lwd=plot.line.thickness)
      box(lwd=plot.line.thickness)
    } else if(pca.visual.flavour == "loadings"){
      biplot(pca.results, 
             col=c("grey70", "black"), 
             cex=c(0.7, 1), xlab="", 
             ylab=PC2_lab, 
             main=paste(graph.title, "\n\n", sep=""), 
             sub=paste(PC1_lab,"\n",graph.subtitle, sep=""),var.axes=FALSE)
    } else if(pca.visual.flavour == "technical"){
      layout(matrix(c(1,2), 2, 2, byrow = TRUE), widths=c(3,1))
      biplot(pca.results, col=c("black", "grey40"), cex=c(1, 0.9), xlab="", ylab=PC2_lab, main=paste(graph.title, "\n\n", sep=""), sub=paste(PC1_lab,"\n",graph.subtitle, sep=""),var.axes=FALSE)
      abline(h=0, v=0, col = "gray60",lty=3)
      # add the subpanel to the right 
      row = mat.or.vec(nc=ncol(pca.results$x),nr=1)
      for (i in 1:ncol(row)){row[,i]<-"grey45"}
      # paint the first two PCS black -- i.e. the ones actually plotted
      row[,1]<-"black"
      row[,2]<-"black"
      barplot(expl.var, col = row, xlab = "Principal components", ylab = "Proportion of variance explained (in %)")
      # set a horizontal dashed line, indicating the psychological 5% barrier  
      abline(h=5, lty=3)
    } else if(pca.visual.flavour == "symbols"){
      # determine labels involved
      labels = c()
      for (c in rownames(pca.results$x)){
        labels = c(labels, gsub("_.*","",c))
      }
      COOR = data.frame(pca.results$x[,1:2], LABEL=labels)
      labels<-c(levels(COOR$LABEL))
      # visualize 
      library(lattice)
      sps <- trellis.par.get("superpose.symbol")
      sps$pch <- 1:length(labels)
      trellis.par.set("superpose.symbol", sps)
      ltheme <- canonical.theme(color = FALSE)      
      lattice.options(default.theme = ltheme)
      pl<-xyplot(data=COOR, x=PC2~PC1, xlab=paste(PC1_lab,"\n",graph.subtitle, sep=""), ylab=PC2_lab, groups=LABEL, sub="", key=list(columns=2, text=list(labels), points=Rows(sps, 1:length(labels))),
             panel=function(x, ...){
             panel.xyplot(x, ...)
             panel.abline(v=0, lty=3)
             panel.abline(h=0, lty=3)
      })
      plot(pl)
    }
  }
}
        

# prepares a list of dendrogram-like structures for a bootstrap consensus tree
# (the final tree will be generated later, outside the main loop of the script)
if (analysis.type == "BCT") {
  mfw.info = paste(mfw.min,"-",mfw.max.original, sep="")
  name.of.the.method = "Bootstrap Consensus Tree"
  short.name.of.the.method = "Consensus" 
  # calculates the dendrogram for current settings
  # 
########################################################################
########################################################################
# compatibility mode: to make one's old experiments reproducible
  if(nj.consensus.tree == TRUE) {
    current.bootstrap.results = nj(as.dist(distance.table))
    } else {
  current.bootstrap.results = as.phylo(hclust(as.dist(distance.table),"ward"))
  }
########################################################################
  # adds the current dendrogram to the list of all dendrograms
  bootstrap.list[[number.of.current.iteration]] = current.bootstrap.results }

  
# establishing the text to appear on the graph (unless "notitle" was chosen)
if(ngram.size > 1) {
  ngram.value = paste(ngram.size,"-grams", sep="")
  } else {
  ngram.value = "" }
  #
  if(titles.on.graphs == TRUE) {
  graph.title = paste(basename(getwd()),"\n",name.of.the.method)
  if(analysis.type == "BCT") {
  graph.subtitle = paste(mfw.info," MF",toupper(analyzed.features)," ",ngram.value," Culled @ ",culling.info,"%\n",
                    pronouns.info," ",distance.name.on.graph," Consensus ",consensus.strength," ",start.at.info, sep="") 
          } else {
          graph.subtitle = paste(mfw.info," MF",toupper(analyzed.features)," ",ngram.value," Culled @ ",culling.info,"%\n",
                    pronouns.info," ",distance.name.on.graph," ",start.at.info, sep="") }
  } else {
  graph.title = ""
  graph.subtitle = "" }


# name of the output file (strictly speaking: basename) for graphs
graph.filename = paste(basename(getwd()),short.name.of.the.method,mfw.info,
                       "MFWs_Culled",culling.info,pronouns.info,
                       distance.name.on.file,"C",consensus.strength,start.at.info, sep="_")
  if(analysis.type == "BCT") {
    graph.filename = paste(basename(getwd()),short.name.of.the.method,mfw.info,
                       "MFWs_Culled",culling.info,pronouns.info,
                       distance.name.on.file,"C",consensus.strength,start.at.info, sep="_") 
  } else {
    graph.filename = paste(basename(getwd()),short.name.of.the.method,mfw.info,
                       "MFWs_Culled",culling.info,pronouns.info, distance.name.on.file,start.at.info, sep="_") 
}

# #################################################
# plotting 
# #################################################

# The core code for the graphic output (if bootstrap consensus tree 
# is specified, the plot will be initiated later)
if(analysis.type != "BCT") {
  if(display.on.screen == TRUE) {
    plot.current.task()
    }
  if(write.pdf.file == TRUE) {
    pdf(file = paste(graph.filename,"%03d",".pdf",sep=""),
            width=plot.custom.width,height=plot.custom.height,
            pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
  if(write.jpg.file == TRUE) {
    jpeg(filename = paste(graph.filename,"%03d",".jpg",sep=""), 
            width=plot.custom.width,height=plot.custom.height,
            unit="in",res=300,pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
  if(write.emf.file == TRUE) {
    if(Sys.info()[["sysname"]] == "Windows") { 
      ### Windows
      win.metafile(filename = paste(graph.filename,"%03d",".emf",sep=""), 
            width=plot.custom.width,height=plot.custom.height,
            res=300,pointsize=plot.font.size)
      plot.current.task()
      dev.off()
      } else {
      ### Linux, Mac
      cat("\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
      cat("EMF file format is not supported by", Sys.info()[["sysname"]],"\n")
      cat("You're suggested to try again with PNG, JPG or PDF.\n")
      cat("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
      }
    }
  if(write.png.file == TRUE) {
    png(filename = paste(graph.filename,"%03d",".png",sep=""), 
            width=plot.custom.width,height=plot.custom.height,
            unit="in",res=300,pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
}
##################################################


# writing distance table(s) to a file (if an appropriate option has been chosen)
if(save.distance.tables == TRUE) {
  distance.table.filename = paste("distance_table_",mfw,"mfw_",current.culling,"c.txt",sep="")
  write.table(file=distance.table.filename, distance.table)
}

# writing the words (or features) actually used in the analysis
if(save.analyzed.features == TRUE) {
  cat(colnames(table.with.all.freqs[,1:mfw]),
     file=paste("features_analyzed_",mfw,"mfw_",current.culling,"c.txt",sep=""),
     sep="\n")
}

# writing the frequency table that was actually used in the analysis
if(save.analyzed.freqs == TRUE) {
  write.table(table.with.all.freqs[,1:mfw],
     file=paste("frequencies_analyzed_",mfw,"mfw_",current.culling,"c.txt",sep=""))
}




}    # <-- the internal loop for(i) returns here
# #################################################

# blank line on the screen
cat("\n")


}    # <-- the main loop for(j) returns here
# ################################################




# bootstrap visualization
if(analysis.type == "BCT") {

# as above, the task to be plotted is saved as a function
if(length(bootstrap.list) <= 2) {
  cat("\n\nSORRY, BUT YOU ARE EXPECTING TOO MUCH...!\n\n",
  "There should be at least 3 iterations to make a consensus tree\n\n")
  } else {
  plot.current.task = function(){ 
        plot(consensus(bootstrap.list, p=consensus.strength),
           type="u",
           font=1,
           lab4ut="axial", 
           tip.color = colors.of.pca.graph)
        title (main = graph.title)
        title (sub = graph.subtitle) }

# The core code for the graphic output... Yes, you are right: you’ve seen
# the same lines above. Instead of blaming us, write better code yourself
# and let us know.
  if(display.on.screen == TRUE) {
    plot.current.task()
    }
  if(write.pdf.file == TRUE) {
    pdf(file = paste(graph.filename,"%03d",".pdf",sep=""),
         width=plot.custom.width,height=plot.custom.height,
         pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
  if(write.jpg.file == TRUE) {
    jpeg(filename = paste(graph.filename,"%03d",".jpg",sep=""),
         width=plot.custom.width,height=plot.custom.height,
         unit="in",res=300,pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
  if(write.emf.file == TRUE) {
    win.metafile(filename=paste(graph.filename,"%03d",".emf",sep=""), 
         width=plot.custom.width,height=plot.custom.height,
         res=300,pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
  if(write.png.file == TRUE) {
    png(filename = paste(graph.filename,"%03d",".png",sep=""), 
         width=plot.custom.width,height=plot.custom.height,
         unit="in",res=300,pointsize=plot.font.size)
    plot.current.task()
    dev.off()
    }
}}


# #################################################
# final cleaning


cat("\n")
cat("removing most of the variables... \n")
cat("type ls() if you want to see what was not removed\n")
cat("if you are going to change the corpus, clean all: rm(list=ls())\n")
cat("\n")


# a list of variables not to be removed
do.not.remove = c("table.with.all.zscores", "table.with.all.freqs",
                  "frequencies.0.culling", "distance.table",
                  variables.not.to.be.removed)

# removing the variables which are not on the above list
list.of.variables = ls()
rm(list=list.of.variables[!(list.of.variables %in% do.not.remove)])


# #################################################


# TO DO:


# Christof: a loop for different start.at values
# Fotis: custom list of files does not work for sample labels

# applicable scenarios:
# 
# 1. MDS, 100, En, pdf, png
# 2. MDS, 1000, 100% culling, En, pdf, png
# 3. PCA, corr., 100
# 4. Cons.Tree
# ...

# common wordlist when the number of full-sized novels >100
# (picking the first 100 by chance? extracting randomly 1M words?,
# extract a number of words, say 50k, from each novel?)

# dendrograms: ward, complete, average, nj
# rooted consensus trees?

# the code for MDS and PCA in different flavors deserves some comments!




Script Guide
------------

Below are the explanations for the purpose and usage of each experimentation script in this folder.

### analysisDynamicsEstimation.m

| Date		| Made for				| Status
=====================================
| ~2018		| single piece analysis	| not up-to-date


### auxFrontiers.m

| Date			| Made for						| Status
=====================================
| late 2018		| perceptual results analysis	| Replaced in paper review

### createfigure_Comparisons.m

| Date			| Made for										| Status
=====================================
| late 2018		| reproducing cosmetic adjustments to graph		| Replaced in paper review


### dumpMidis.m

| Date				| Made for                                  | Status
=====================================
| feb 2019          | Quickly producing MIDIs for setup tests	| could be reworked


### exaggerated_notes_arff.m

| Date			| Made for					| Status
=====================================
| early 2019	| create note-level ARFFs	| *up-to-date*


### expertDBLeaveOneOut.m

| Date		| Made for				| Status
=====================================
| ~2017		| early model analysis	| legacy


### generateDynViolinRT.m

| Date			| Made for						| Status
=====================================
| late 2017		| SkyNote reference creation	| not up-to-date


### GenPercepTestMid.m

| Date			| Made for						| Status
=====================================
| late 2018		| perceptual test preparation	| up-to-date (needs test)

This script generates the midi files for blind testing and also produces
the boxplot of leave-one-piece-out mean absolute error per phrase.

### learningCurve.m

| Date				| Made for                      | Status
=====================================
| feb 2019          | extra plot in Frontiers paper	| *up-to-date*

This script runs a cross-validation analysis with an incrementally 
larger dataset built from randomizing all phrases read from the files in
the folder given by inputFolder. In the end, a plot of the evolution of 
mean absolute error as the dataset grows is supposed to provide an idea on
weather the amount of data provided is enough.

### main.m

| Date				| Made for				| Status
=====================================
| early 2017		| early model analysis	| legacy


### printfigpdf.m

| Date		| Made for				| Status
=====================================
| 2017?		| automatic PDF figure	| needs review


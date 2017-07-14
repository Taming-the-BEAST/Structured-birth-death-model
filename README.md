---
author: Denise Kühnert, Jūlija Pečerska
level: Professional
title: Structured birth-death model
subtitle: Population structure using the multi-type birth-death model
---


# Introduction

In this tutorial we will use the [BEAST2](http://www.beast2.org/)
[bdmm](https://github.com/denisekuehnert/bdmm) package to perform a Bayesian
phylogenetic analysis of an influenza data set using the multi-type birth-death
model.

<!--(Note that both the structured coalescent and the multi-type birth-death model are tree priors implemented in BEAST2. Both of them utilize the multi-type tree structure of the MultiTypeTree package. While the structured coalescent is part of the MultiTypeTree package, the multi-type birth-death model has its own package bdmm (aka birth-death migration model).)-->


The data set used in this tutorial is a thinned 60 sequence subset of the
980 sequence H3N2 influenza data set used in the publication {% cite Vaughan2014 --file Structured-birth-death-model/refs.bib %}, which in turn was
assembled from publicly-available data sets provided by various authors on
[GenBank](http://www.ncbi.nlm.nih.gov/genbank/).

## Software Requirements

### BEAST2 - Bayesian Evolutionary Analysis Sampling Trees 2

This tutorial uses the [BEAST2](http://www.beast2.org/) version 2.4.7. The software package also includes BEAUTi, the graphical BEAST2 interface which we will use to set up the analysis.

### Tracer

To summarise the inference results you will also need a recent version of [Tracer](http://tree.bio.ed.ac.uk/software/tracer/).

<!-- and an up-to-date version of 
[Google Chrome](http://www.google.com/chrome) or
[Mozilla Firefox](https://www.mozilla.org/en-US/firefox/).-->


## Installing the `bdmm` package

You can easily install the `bdmm` package via BEAUti's package manager.  To do this, follow these steps:

1. Start BEAUti;
2. In the application menu, `File > Manage packages`.
3. Find `bdmm` in the list of packages shown, select it and then click `Install/Upgrade`.

The BEAUTi window should look similar to what is shown in [Figure 1](#fig:install-bdmm).
Note the actual version of `bdmm` may differ from the version shown in the figure, which is perfectly normal.
Also note that `bdmm` depends on the `MultiTypeTree` package. BEAUTi will install the dependencies automatically once you select to install `bdmm`.

<figure>
	<a id="fig:install-bdmm"></a>
	<img style="width:75%;margin:auto;display:block;" src="figures/1-install-bdmm.png" alt="">
	<figcaption style="width:75%;margin:auto;display:block;">Figure 1: Install bdmm.</figcaption>
</figure>
<br>

Finally, **restart BEAUti.**  The restart is necessary for the packages to be successfully installed.

# Setting up the analysis using BEAUti

## Loading the Template

A BEAUTi template defines the basic structure and contents of your XML configuration file.
By default BEAUTi will construct an XML file with standard uncoloured BEAST trees, however `bdmm` uses coloured trees which are defined in the `MultiTypeTree` package.
To use the appropriate template for the configuration file, select `File > Template > MultiTypeBirthDeath`, as shown in [Figure 2](#fig:choose-bdmm).

<figure>
	<a id="fig:choose-bdmm"></a>
	<img style="width:100%;" src="figures/2-choose-bdmm-template.png" alt="">
	<figcaption>Figure 2: Load the MultiTypeBirthDeath template.</figcaption>
</figure>
<br>


## Loading the data

Once the template is loaded, we can load in our example sequence data.  In our case, this data is stored in a FASTA file, the first few lines of which look like this (the sequences have been truncated for better readability):

```
> EU856841_HongKong_2005.34246575
-----------GGGATAATTCTATTAACCATGAAGACTATCATTGCTTTGAGCTACATTT...
> EU856989_HongKong_2002.58356164
--CAAAAGCAGGGGATAATTCTATTAACCATGAAGACTATCATTGCTTTGAGCTACATTT...
> CY039495_HongKong_2004.5890411
------------------TTCTATTAACCATGAAGACTATCATTGCTTTGAGCTACATTC...
> EU856853_HongKong_2001.17808219
---------------------TATTAACCATGAAGACTATCATTGCTTTGAGCTACATTC...
> CY010084_NewZealand_2005.62739726
---------------------TATTAACCATGAAGACTATCATTGCTTTGAGCTACATTC...
> CY007387_NewZealand_2004.63287671
---------------------TATTAACCATGAAGACTATCATTGCTTTGAGCTACATTC...
> CY012432_NewZealand_2000.81643836
---------------------------CCATGAAGACTATCATTGCTTTGAGCTACATTT...
```

The lines beginning with ">" are labels for the sequences immediately
following.  In general, these labels have no special format, but in this file
each label is an underscore-delimited triple.  The first element of each triple
is the GenBank accession number of the sequence, the second is the geographical
region from which it was sampled, and the third is the time at which it was
sampled measured in calendar years or fractions thereof. 

In this tutorial we will be using the influenza sequence data which can be found in the `examples` folder of the `MultiTypeTree` package. 
To make it easier to find when loading the alignment, you can optionally set the working directory of BEAST2 to `MultiTypeTree`.
This will make BEAUTi open the appropriate package folder when you look for the alignment.
To set the working directory, select `File > Set working dir > MultiTypeTree`, as shown in [Figure 3](#fig:working-dir).

<figure>
	<a id="fig:working-dir"></a>
	<img style="width:100%;" src="figures/3-set-working-dir.png" alt="">
	<figcaption>Figure 3: Optional step: set the working directory to MultiTypeTree.</figcaption>
</figure>
<br>

To load the file, select `File > Add alignment`.

This will open a file selection dialog box.  The example influenza sequence data
file is named `h3n2_2deme.fna`.
Assuming you have followed the previous step to set the working directory, this can be found in the `examples/` directory shown when the file selection dialog box appears.
In case you have not followed the previous step you will have to locate the folder containing the `MultiTypeTree` package and look for the `examples/` folder there.

Once the sequence file is loaded, your BEAUti screen should look similar to what is shown in [Figure 4](#fig:alignment).

<figure>
	<a id="fig:alignment"></a>
	<img src="figures/4-alignment-loaded.png" alt="">
	<figcaption>Figure 4: The alignment loaded into BEAUti.</figcaption>
</figure>
<br>

## Setting up dates

Once the data is loaded, the next step is to specify the times at which the sequences were sampled:

1. Select the `Tip Dates` panel.
2. Check the `Use tip dates` checkbox.
3. Click the `Auto-configure` button at the top-right of the panel.
This opens a dialog that allows sample times to be loaded from a file or inferred (guessed) from the sequence labels.
4. Because the times are included as the last element of the underscore-delimited sequence names, choose the `use everything` radio button and select `after last` from the drop-down menu. The default delimiter is already the underscore, so there is no need to change that.

The date parsing setup will look as shown in [Figure 5](#fig:tip-dates).

<figure>
	<a id="fig:tip-dates"></a>
	<img style="width:75%;margin:auto;display:block;" src="figures/5-tip-dates.png" alt="">
	<figcaption style="width:75%;margin:auto;display:block;">Figure 5: Guessing the tip dates.</figcaption>
</figure>
<br>

After clicking `OK` you should find that the tip date table is filled with
times that match those in the sequence headers, and that the last column of the
table contains heights, i.e. times before most recent sample, calculated from the times.
The BEAUTi panel should look as shown in [Figure 6](#fig:tip-dates).

<figure>
	<a id="fig:tip-dates-set"></a>
	<img src="figures/6-tip-dates-set.png" alt="">
	<figcaption>Figure 6: Sampling dates as seen in BEAUti.</figcaption>
</figure>
<br>

## Setting up locations

Now that we've specified the sampling times, we move on to specifying the sampling locations.
To do this, we follow a very similar set of steps to those we used to set the sample times:

1. Select the `Tip Locations` panel. You'll find that the locations are already filled with a single default value – `NOT_SET`. 
2. Click the `Guess` button at the top-right of the panel. This opens the same dialog that we saw in the previous section when setting up the dates.
3. The locations are included as the second element of the underscore-delimited sequence names.
Therefore we choose the `split on character` radio button and select group `2` from the drop-down menu.
Note again that the underscore character is already chosen as the delimiter.

The location parsing setup will look as shown in [Figure 7](#fig:tip-types).

<figure>
	<a id="fig:tip-types"></a>
	<img style="width:75%;margin:auto;display:block;" src="figures/7-tip-types.png" alt="">
	<figcaption style="width:75%;margin:auto;display:block;">Figure 7: Guessing the locations.</figcaption>
</figure>
<br>

After clicking `OK` you should find that the tip location table is filled with locations that match those in the sequence titles.
The BEAUTi panel should look as shown in [Figure 8](#fig:tip-types-set).

<figure>
	<a id="fig:tip-types-set"></a>
	<img src="figures/8-tip-types-set.png" alt="">
	<figcaption>Figure 8: The locations in BEAUti.</figcaption>
</figure>
<br>

## Setting the substitution model

For this analysis, we will use the HKY substitution model with 4 gamma categories and estimated base frequencies.
To configure this in BEAUti, switch to the `Site Model` panel.
First, we need to set up the rate category count.
To approximate the continuous gamma rate distribution BEAST2 uses the discrete gamma distribution, where sites are divided into k equally probable rate categories.
In general, 4-6 categories work well for most datasets, while having more categories involve a lot of computation at little precision gain, so we set the `Gamma category count` to 4.
We would also like to estimate the `Shape` parameter, which describes the chape of the continuous gamma distribution we approximate.
To do so, we need to set it to a non-zero value (e.g. the default 1.0) and tick the `estimate` chackbox.
While the gamma categories account for rate variation, allowing some sites to have an evolutionary rate of 0 can imporove fit to real data.
Thus we we need to set the `Proportion Invariant` parameter to a non-zero value (e.g. 0.5) and tick the `estimate` chackbox.

Next, to set up the substitution model, select `HKY` from the drop-down menu (the default option is `JC69`).
We would like to estimate the kappa parameter of HKY, so we leave the `Kappa` at the default value of 2.0 and leave the `estimate` checkbox checked.
As we would also like to estimate nucleotide frequencies, so we leave the `Frequencies` parameter at the default value (`Estimated`).
The BEAUti panel should now look as shown in [Figure 9](#fig:site-model).

<figure>
	<a id="fig:site-model"></a>
	<img style="width:100%;" src="figures/9-sitemodel.png" alt="">
	<figcaption>Figure 9: Setup of the site model.</figcaption>
</figure>
<br>

Note that the `Substitution rate` defined on this panel should not be estimated - we use the `Clock rate` defined in the `Clock Model` panel to
determine the average per unit time rate of sequence evolution.
This way, the `Substitution rate` is not actually a rate, but rather a rate multiplier that we fix to 1 for better parameter identifiability.

## Setting the clock model

To speed up the analysis we will assume a strict clock for this small dataset.
Since our data comes from a single epidemic it is not a far-fetched assumption, however to make sure this assumption holds model comparison is necessary.
We will not do this here, but selection of clock model for a different, real analysis should not be taken lightly.
Since our alignment contains sequences sampled at different times and those times are measured in years, we must use a real clock rate expressed in units of expected substitutions per site per year.
Usually the precise value is unknown and so the default behaviour of BEAUti is to assume this rate has to be estimated.
To speed up mixing we set the starting value of the `Clock rate` to 0.005, which we know from research to be much closer to the truth than the default value of 1.
The `Clock Model` panel should now look as shown in [Figure 10](#fig:strict-clock).

<figure>
	<a id="fig:strict-clock"></a>
	<img src="figures/10-strict-clock.png" alt="">
	<figcaption>Figure 10: Fix the clock rate to speed up mixing.</figcaption>
</figure>
<br>

## Adjusting Priors

### Setting up the `bdmm` tree prior

`bdmm` defines a prior on the multi-type tree distribution.
Thus it is particularly important for the analysis that we properly set up the priors.
First, let's talk about the values that need to be set on the `Priors` panel.
The first panel that you see at the top is the tree prior.

`bdmm` is a model that can be used to explain data that is clearly divided into separate partitions, or demes.
The demes can be geographical locations, as in our example, but the sequences can also be separated through other means than that, e.g. by a specific drug resistance mutation (strains can develop/lose drug resistance and thus move between demes, but can not transfer between demes otherwise), or location in the body (for example, for localised infections caused by the same agent).
In this dataset we have strains from 2 different locations, New Zealand and Hong Kong, so the `Number of demes` should be set to 2, which luckily also is the default value.
Next, `bdmm` lets you estimate the `Reproduction number per type` and the `BecomeUninfectiousRate per type`.
This will let us see the differences in reproduction fitness and speed of recovery between the two locations, so we leave the `estimate` checkboxes checked.
We can leave the starting values at default as it will not influence the inference a lot.

The next important thing one should take care of is setting the sampling proportions appropriately.
In general, the trees that we build go back in time much further than the first sample that we have.
If we set the same sampling proportion for the whole time period from our estimated tree origin to the time of the last sample, we will most likely run into trouble, as `bdmm` will try to produce a tree that has the same sampling proportion for the whole time, but no samples in the past and a lot of samples towards the present.
In order to remove that bias from the trees, we need to make sure that we only have non-zero sampling starting from the first sample date.
To do so, let's look at the `SamplingProportion per type` field.
You will see that it has 4 values, which correspond to two values per type, lets call them [v1,v2,v3,v4].
v1 and v2 are the values for the first and second time interval for the first deme, and v3 and v4 are the values for the second deme.
Thus, to do what we want we need to set the values v1 and v3 to zero.
Because BEAST2 will use scalers to sample new values for the sampling proportions, the values which we set to 0 will remain so.
Next, we also need to set the `Sampling change time` to the time slightly before the first sample.
If we look back at the `Tip dates` panel, we can see that our oldest sample is the one labelled as `EU856904_HongKong_2000.09863014`, for which the height, or the length of time from the first sample and the last, is 5.569863.
We set the sampling change time in time units from the most recent sample and we need to make sure we include the first sample, thus we set the `Sampling change time` to 5.57, which is the height of the first sample rounded slightly up.
The final setup of the tree prior can be seen in [Figure 11](#fig:tree-prior).

<figure>
	<a id="fig:tree-prior"></a>
	<img src="figures/11-tree-prior.png" alt="">
	<figcaption>Figure 11: Set the change time for the sampling proportion so it is zero before the time of the first sample.</figcaption>
</figure>
<br>

<!--When you expand the tree prior element, you can change the condition on survival setting. We'll leave the box checked. 

<figure>
	<a id="fig:"></a>
	<img src="figures/9b-condition.png" alt="">
	<figcaption>Figure 11: Condition on survival.</figcaption>
</figure>
<br>-->

#### What if you have more demes?

First things first, for an analysis with more demes you need to set the `Number of demes` to the appropriate value, e.g. N, that actually corresponds to the number of demes in the dataset.
When you do that, the dimensions of the parameters `Reproduction number`, `BecomeUninfectiousRate`, `SamplingProportion` and `Migration rates` will change.
The `Reproduction number` and the `BecomeUninfectiousRate` will have as many values as you have demes.
The dimensionality of the `SamplingProportion` will be the number of demes times 2, so in case you have 4 demes your sampling proportion will need 8 values.
You can view this parameter as a matrix of 2 x N values, which is flattened by row.
The first column reflects the sampling rate before first sample and all of the values in it should be set to 0.
The `Sampling change time` obviously does not change dimensionality, but has to be set to the appropriate time for your dataset.
Finally, the `Migration rates` will have N * (N - 1) entries.
As one can imagine, the matrix should have the dimensions of N * N, however since there is no migration from a deme to itself (values on the diagonal of the matrix), we subtract N values from the dimensionality, getting N * (N - 1).

### Setting up other priors

By default, BEAST2 provides you with a prior distribution for each of the parameters of your model.
This is done because otherwise BEAUTi will have a hard time displating all of the parameters without any settings provided.
Unfortunately, this means that some priors are very generic, and, moreover, some priors are in fact, improper – the distribution does not integrate to one.
This means that while the default setup might work and the runs will eventually mix, it can happen that the values are 

So, let us go through the important parameters and set priors according to the information we have about our dataset.
The first important parameter is R0.
In epidemiology, the basic reproduction number, R0, of an infection is the number of secondary cases one case generates on average over the course of its infectious period, in an otherwise uninfected population.
Even though we do not have any information on R0 in the particular outbreak, infections rarely have R0 > 10, so we can set an upper limit on the sampled values.
To do so, in the line denoted `R0.t:h3n2_2deme` click the button captioned `initial = [2.0, 2.0] [0.0, Infinity]` to get a pop up settings window (see [Figure 12](#fig:R0-prior)), where you can set the upper value.
Other than that, the current prior sets ther median value of the ditribution to e<sup>0</sup> = 1, which will fit the endemic case of influenza.
<!--Should we also have multiple dimensions for R0 here?-->

Next, we should adjust the prior for the rate of clearing the infection, which is labelled as `becomeUninfectiousRate.t:h3n2_2deme`.
The value of the rate, say x, is the reciprocal of the average time a person with influenza is infectious, 1/x.
Given that we know little about this particular epidemic, let us just assume that a person is infectious for at least 1 day and set the prior for this parameter to `Uniform` with the upper value of 365 (see Figure 13](#fig:bUR-prior)).

For the purpose of this tutorial and given that we know little about the outbreak in question, we will leave the other priors on the default values, but feel free to through the other priors yourself and verify their sensibility.

<figure>
	<a id="fig:R0-prior"></a>
	<img src="figures/12-R0-prior.png" alt="">
	<figcaption>Figure 12: Set the prior for the R0.</figcaption>
</figure>
<br>


<figure>
	<a id="fig:bUR-prior"></a>
	<img src="figures/13-bUR-prior.png" alt="">
	<figcaption>Figure 13: Set the prior for the rate of recovery.</figcaption>
</figure>
<br>

We will use the default set-up for the MCMC and save our file as usual.

## Saving the configuration

Once you are done wuth setting all the appropriate parameters, you can save the configuration file.
For now we will leave the `MCMC` panel parameters as they are by default.

# Running the analysis using BEAST

To run the analysis, simply start BEAST 2 in the manner appropriate for your platform, then select the configuration file you generated in the last section as the input.
Unfortunately, this particular run will take quite some time to mix, e.g. on a 2.5 GHz i7 MacBook Pro it takes about 6 hours for 10'000'000 samples.
Feel free to run it and observe the results, but for the purpose of finishing the tutorial in a reasonable time, check out the provided log file to see the results.

# Analyzing the results

The results of the analysis primarily consist of two parts:

1. The parameter log, which is written to the file `h3n2-bdmm.log`.
2. The tree log, which is written to `h3n2-bdmm.h3n2_2deme.trees`.

In addition, the file `h3n2-bdmm.h3n2_2deme.map.trees` contains the running
estimate of the MAP tree as a function of MCMC step number, while the file
`h3n2-bdmm.h3n2_2deme.typedNode.trees` is the TreeAnnotator-compatible file
we'll use to assemble a summary tree.

## Parameter log file analysis

We can use the program [Tracer](http://tree.bio.ed.ac.uk/software/tracer/) to view the parameter log file. To do this, start Tracer and then press the `+` button in the top-left hand corner of the window (under `Trace files`).
Select
the log file for this analysis (`h3n2_2deme.log`) from the file selection dialog box.
The `Traces` table will then be populated with parameters and summary
statistics corresponding to our multitype birth-death analysis.

Important traces are:
* `R0.t:h3n2_2deme1` and `R0.t:h3n2_2deme2`: These give the effective reproduction numbers for deme 1 (Hongkong) and 2 (New Zealand), respectively.

* `rateMatrix.t:h3n2_2deme1`	`rateMatrix.t:h3n2_2deme2`: These give the (per lineage per year) migration rates from deme 1 to 2 and vice versa.

* `Tree.t:h3n2_2deme.count_HongKong_to_NewZealand`: these give the actual number of ancestral
  migrations from HongKong to New Zealand.

The panels tabs at the top-right of the window can be used to display one or
more selected traces in various ways.  For example, selecting the two R0 traces and choosing the "Marginal prob distribution" panel results in the following useful comparison between the sampled population size
marginal posterior distributions:

<figure>
	<a id="fig:"></a>
	<img src="figures/tracer-R0.png" alt="">
	<figcaption>Figure 13: Estimated {% eqinline R_0 %} marginal posteriors.</figcaption>
</figure>
<br>

Note that some of the ESS values are still less than 200 - the arbitrary
threshold for acceptability. If this analysis were part of a serious study you
would want to run the chain for another few million iterations to improve
these. (In BEAST 2, analyses can be resumed - the samples you've already
acquired will not be wasted.) For the purposes of this tutorial, however, these
values are acceptable.

## Tree log visualization

The popular phylogenetic tree visualizer
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/) can be used to visualize
the sampled trees.  Be warned, however, that FigTree currently
takes an extremely long time to load even relatively small (a few megabyte)
MultiTypeTree logs.

For this reasons we suggest using [IcyTree](http://tgvaughan.github.io/icytree)
to view tree log files and maybe switching to FigTree to visualize summary
trees as discussed in the next section.  (Also, IcyTree can be used to export
individual trees from a large log file for subsequent viewing using FigTree.)
IcyTree is a tree viewer that runs in a web browser.  It runs best under recent
versions of [Google Chrome](http://www.google.com/chrome) and [Mozilla
Firefox](https://www.mozilla.org/en-US/firefox/) (in that order).

To view MultiTypeTree log files using IcyTree, simply navigate to the IcyTree
web page, select "Load from file" from the "File" menu, then select one of the
tree log files using the file selection dialog. Once the file is loaded you
will see the first tree it contains.  In order to select a different tree, move
the mouse pointer over the box in the lower-left corner of the window.  This
box will expand to a small dialog containing buttons allowing you to navigate
between trees. The '<' and '>' buttons move in steps of 1 tree, while '<<' and
'>>' move 10% of the tree file per click.  You can also directly enter the
index of a tree.  (Note that there are keyboard shortcuts for almost all
commands in IcyTree and that these can be found by selecting "Keyboard
shortcuts" from the "Help" menu.)

Initially the trees edges will be uncoloured.  To colour the edges according to
the edge type, open the "Style" menu, navigate to the "Colour edges by" submenu
and select "type". A legend and axis can be added by choosing "Display legend"
and "Axis > Age" from the same menu.

The following shows the final tree of `h3n2-bdmm-v2-samplingPrior.h3n2_2deme.map.trees` in IcyTree, which
represents our sampled estimate of the MAP multi-type tree:

<figure>
	<a id="fig:"></a>
	<img src="figures/icyTreeMAP.png" alt="">
	<figcaption>Figure 14: The MAP multi-type tree in IcyTree.</figcaption>
</figure>
<br>

While IcyTree is useful for rapidly visualizing the results of an analysis, it
is not nearly as feature-rich as FigTree and not as capable for producing
publication-quality graphics.  Happily, however, IcyTree can extract single
trees from larger log files. Simply navigate to the desired tree, open the
"File" menu, choose the "Export tree as" submenu and select "NEXUS file". (It
is important to select "NEXUS" instead of "Newick" as the Newick format does
not support the annotations that MultiTypeTree uses to mark the edge types.)

## Producing a summary tree using TreeAnnotator

While it is tempting to view the MAP tree shown above as the primary result of
the phylogenetic side of our analysis it is very important to remember that
this is only a point estimate and says nothing about the uncertainty present in
the result.  This is an important drawback, as we have done a full Bayesian
analysis and have access to a large number of samples from the full posterior
in the tree log files. The MAP tree discards almost all of this information.

We can make better use of our raw analysis results by using the TreeAnnotator
program which is distributed with BEAST to analyze the
`typedNode` trees which were produced by our MCMC run.  To do this,
simply load TreeAnnotator and select the `typedNode` tree file as
the input file and `h3n2-bdmm-v2-samplingPrior.h3n2_2deme.summary.trees` as the output file.  Select "Mean
heights" from the "Node heights" menu and set the burn-in percentage to 10:

<figure>
	<a id="fig:"></a>
	<img style="width:50%;" src="figures/treeannotator.png" alt="">
	<figcaption>Figure 15: Use TreeAnnotator to produce a summary tree.</figcaption>
</figure>
<br>

Pressing the "Run" button will now produce an annotated summary tree.

To visualize this tree, open IcyTree once more (maybe open it in a new browser
tab), choose File->Open, then select the file
`h3n2_2deme.h3n2_2deme.summary.tree` using the file selection dialog. Follow
the instructions provided for the MAP tree above to colour the tree by the
"type" attribute and add the legend and time axis. In addition, open the Style
menu and from the "Node height error bars" sub-menu select "height_95%_HPD" to
add error bars to the internal node heights. Also, open the Style menu and from
the "Edge opacity" sub-menu select "type.prob". This will cause the edges to
become increasingly transparent as the posterior probability for the displayed
colour decreases.

Once these style preferences have been set, you should see something similar to
the following:

<figure>
	<a id="fig:"></a>
	<img src="figures/icyTreeSummary.png" alt="">
	<figcaption>Figure 16: The summary tree in IcyTree.</figcaption>
</figure>
<br>

Here we have a full consensus tree annotated by the locations at coalescence
nodes and showing node height uncertainty, with the widths of the edges
representing how certain we can be of the location estimate at each point on
the tree. This is a much more comprehensive summary of the phylogenetic side of
our analysis.

One thing to pay attention to here is that the most probable root location is
given by the summary tree to be Hong Kong (under our model which assumes that
only Hong Kong and New Zealand exist). By hovering the mouse cursor over the
tiny edge above the root will bring up a table in which posterior probability
of the displayed root location (`type.prob`) can be seen to be approximately
90%. The analysis therefore strongly supports a Hong Kong origin over a New
Zealand origin.  

<!--[Very useful final notes from Tim](https://github.com/CompEvol/MultiTypeTree/wiki/Beginner%27s-Tutorial-%28short-version%29#final-notes)-->


----

# Useful Links

- [Bayesian Evolutionary Analysis with BEAST 2](http://www.beast2.org/book.html) {% cite BEAST2book2014 --file Structured-birth-death-model/refs.bib %}
- [Multi-type birth-death process package](https://github.com/denisekuehnert/bdmm) {% cite Kuhnert2016 --file Structured-birth-death-model/refs.bib %}
- BEAST 2 website and documentation: [http://www.beast2.org/](http://www.beast2.org/)

----

The content of this tutorial is based on the [Structured Coalescent tutorial](https://github.com/CompEvol/MultiTypeTree/wiki/Beginner's-Tutorial-(short-version)) by Tim Vaughan.

-----

# Relevant References

{% bibliography --cited --file Structured-birth-death-model/refs %}
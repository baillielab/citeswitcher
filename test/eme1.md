--- 
title: Mechanistic sub-study within the HTA-funded A2B trial testing immunomodulatory mechanisms of alpha_2-andrenoreceptor antagonists. 
date: \today 
bibliography: eme1.bib
...

[comment]: <> (NOTES: 
Baillie JK, Walsh TS, Campbell L, Shankar-Hari M, MacLullich A, Singer M, Gordon A, McAuley D for A2B investigators
-----
REQUIREMENTS:
font size = 11 Arial
20 pages max
header with project ref number
footer with page numbers
pdf
References seperately 
-----
TO DO 
Correct immunomodulatory and organ dysfunction throughout.
)

[comment]: <> (Please use this format: [PMID:xxxx] for citations)

# Background 

In mid-2018 we will start an HTA-commissioned (16/93/01) three-arm, parallel-group, open-label, randomized clinical trial to establish whether either of the alpha_2-adrenoreceptor agonists (alpha_2-agonists) dexmedetomidine or clonidine are superior to usual care (propofol) as primary sedative agents in mechanically ventilated patients (‘A2B trial’). The question posed by the HTA commission was would alpha_2-agonists compared to propofol improve clinical outcomes due to their potentially superior sedative properties, albeit with uncertain safety and cost-effectiveness. 

However, alpha_2-agonists are also known to have extensive immune-modulating effects in both *in vitro* and *in vivo* model systems relating to critical illness[PMID:24856796; PMID:19691839]. Since immune mechanisms precipitate and perpetuate organ failure in critical illness, immunodulatory effects may mediate the potential beneficial clinical effects of these agents in critically ill patients.

If the A2B trial shows any clinical effect from treatment with alpha_2-agonists, this important mechanistic question will immediately follow: is the effect due to the sedative or immunodulatory effects of these drugs? Using a single RNA-derived composite measure of alpha_2-agonist immunodulation, we can provide robust statistical evidence to answer this question.

The randomisation in the A2B trial provides a unique instrumental variable to evaluate the causal effects of inflammatory pathways on key outcomes in human critical illness.

![Directed acyclic graph demonstrating proposed causal infrence. Black circles show exposures, white circles show outcomes. Solid black arrows indicate known causal relationships.](img/a2b_schematic.png)

## Modifiable immune activation in critical illness 

Organ failure is the defining feature of critical illness. Immune mechanisms are believed to play a causative role in organ failure from multiple distinct precipitants, including sepsis, burns, trauma, haemorrhage, pancreatitis and others. In each of these conditions, inflammatory organ damage occurs at sites *remote* from the initial injury. Our recent systematic review [ref: systematic review] demonstrates that alpha_2-agonists have direct effects on both quiescent and activated immune cells, and induce complex modulation of inflammatory responses in animal models. 

### Quantifying alpha_2-agonist-mediated immunomodulation

In order to demonstrate causal relationships between immune activity and clinical outcomes, it will be necessary to quantify the relevant component of immune activation. Specifically, we need to quantify the component of immune activity that is modifiable by alpha_2 agonists. 

For the following two reasons, this problem requires a novel approach. Firstly, immune activation in critical illness is a complex multidimensional phenotype with a high degree of inter-individual variation[PMID: 24855243]. Secondly, alpha_2-agonists have cell type-specific modulatory effects on both innate and adaptive immunity *in vivo*[PMID:24856796; PMID:19691839]. In order to draw a robust causal inference with optimal statistical power we will decompose these complex signals into a single primary composite measure: alpha_2-agonist-mediated immunomodulation. The combined cost of measuring cytokines, other mediators, and separating immune cell types to measure gene expression, would be prohibitive in a study on the scale proposed here.

To overcome this problem we have:
1. completed a systematic review of *in vitro* and *in vivo* studies of alpha_2-agonist effects on immune function. From this review, we have identified speficic gene expression signatures defining the alpha_2-agonist-modifiable components of immune activation [figure; see below].
2. initiated a collaboration with Prof. Tom Freeman, Roslin Institute, who has recently completed a comprehensive evaluation of RNA-based immune cell signatures (ImSigDB - the immune signatures database), identifying key covariates to enable systematic deconvolution of cell type-specific signals from whole blood.
3. Devised a novel robotics-driven protocol to enable high-throughput cap analysis of gene expression (CAGE) sequencing at low cost, enabling the identification of highly cell type-specific RNA transcriptional signatures.

Together, these steps enable the identification of a single, RNA-derived signal defining the alpha_2-agnonist-mediated immunomodulation in whole blood from critically-ill patients. This can be obtained at low cost during treatment with an alpha_2-agonist or placebo.

[what if alpha_2-agonist effect is mediated through immune cell number and we eliminate this signal during deconvolution? ==> answer, we'll do a raw analysis too.]

### Sedation
[Tim - could you write a bit about quantifying sedation here?]

## Outcomes

In the A2B trial our primary outcome is the clinically and economically important ‘time to successful extubation’, which is strongly affected by underlying inflammatory disease processes [ref]. Inflammatory pathways also have important effects on key secondary outcomes in the trial, including delirium[PMID:22884900], mortality, and cognitive decline[PMID:28114436].

### Time to extubation

### Delirium

Delirium is a complex multifactorial clinical phenotype, with high prevalence during critical illness, and is a key secondary trial outcome in A2B. Delirium is strongly associated with multiple clinically important adverse outcomes from critical illness, including mortality and subsequent cognitive impairment. Significant evidence supports a lower delirium incidence with alpha_2-agonist based sedation. Inflammation is strongly associated with delirium, but intervention studies in humans - necessary to infer causation - are lacking. Elucidating the mechanism by which alpha_2-agonists alter delirium could give novel insights into delirium pathogenesis and potential future therapies.

[ref Gerrard delirium endotypes]

# Hypotheses 

1. We hypothesise that clinical outcomes are significantly mediated through modification of systemic inflammation. 

[comment]: <> (specific eQTLs specific CAGE tags)


# Aims 

Test in the A2B trial whether inflammatory signals are causally related to a reduction in time to extubation, delirium, or mortality.

## Study design 

### Population
We will select a subset of patients recruited to the A2B trial in whom there is evidence of systemic inflammation and a high severity of illness. 
[+++ details]

### Exposures

#### Definition of composite measure of alpha_2-agonist immunomodulation

In order to inform a hypothesis-testing study in patients, we have collated published studies of effects of alpha_2-agonists on transcript and cytokine production in immune cells from any species and systematically collapsed gene signals onto human gene names using publicly-available annotation and orthology data, as in our previous work [PMID:22451944]. We employed a novel  crossvalidation algorithm (meta-analysis by information content, MAIC) to systematically combine and evaluate data form different sources. 

### Sample acquisition

In order to detect the biological effect of treatment we will obtain a single blood sample at 48-72h after randomisation from all consenting patients in the A2B trial. In order to cost-effectively measure multiple transcripts arising from distinct genomic regions in different cell types, CAGE RNA sequencing will be performed in a randomly-selected subset of 300 patients, 100 from each of the study groups. 

### Outcomes

Time to extubation
Delirium
Mortality


## Causal inference

For each outcome, we will quantify evidence for causality of both sedative and immunomodulatory effects (Figure 1) by combining (a) statistical evidence for a difference in each of these exposures between treatment and control groups and (b) statistical evidence for an association with outcome after adjusting for other measured characteristics. 

After correction for multiple comparisons we will identify causative evidence for specific outcomes. Inflammatory signals identified in this instrumental variable analysis will either lie on the causal pathway between alpha_2-agonists and key outcomes, or will be biomarkers (in this biological context) for causative processes occurring outside of the scope of our measurements (e.g. unmeasured mediators in blood, immune cells in solid organs, direct effects on the central nervous system).

### Serial samples vs. single sample

Additional information relevant to our hypothesis test could be obtained by obtaining serial samples, for example to evaluate the change from baseline after initiation of treatment with alpha_2-agonists. For a finite study budget, this would decrease by half the number of patients that we can study. Our power modelling [see below] indicates that the anticipated gain in signal:noise ratio would not overcome the reduction in study power. We have therefore concluded that a single sample, 48 hours after initiation of treatment, will provide the greatest probability of detecting an effect.

# Deliverables

1. Determine whether clinical differences between treatment groups are mediated by anti-inflammatory, rather than sedative, effects of alpha_2-agonists

# References

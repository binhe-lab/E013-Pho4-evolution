---
title: 
author: 
date: 
---

# 2022-11-03, Yuning, probe design, Lindsey reply adding Cg probes

Note that when comparing the result data to Yuning's letter, most categories match, except for the "Sc genome Pho4/Cbf1 binding sites". I found a total of 2904 probes with the name "sc_Probes_Cbf1_Pho4_uniq_mutMax". This fits my expectation as there are a total of ~800 CACGTG motifs in the yeast genome. After some data transformation, I identified 68 consensus motif probes, roughly corresponding to the number of ChIP peaks for ScPho4, and also 417 "nonconsensus" probes. I'm not sure how those were chosen. Need to look into Yuning's paper.

Yuning:

> Hi Dr. He,
> 
> I’ve let Raluca know. I think Raluca is doing a final check on the array design (we usually check several rounds before finally placing the order). So at the same time, @Lindsey can you send us a simple spreadsheet for the design with information as follows? Here I already listed the information for probes from yeast genome. So if you can add the info for the Cg probes that you have, that’d be perfect.
> 
> Probe set 1: MITOMI control probes
> Number of probes: 4000x6(replicates)=24000
> Probe name prefix: MitomiProbe_mut
> 
> Probe set 2: Sc genome Pho4/Cbf1 binding sites
> Number of probes:
> 5711x6=34266
> Probe name prefix: sc_Probes_Cbf1_Pho4_uniq_mutMax
> 
> Probe set 3: Sc genome control probes (accessible region but not bound by Pho4)
> Number of probes: 150x6=900
> Probe name prefix: sc_NegCtrl
> 
> Best,
> Yuning

Lindsey:

> Hi Yuning,
> 
> I put the information for the probes below. Let me know if you need anything else.
> 
> Probe set 4: Cg genome Pho4 consensus binding sites
> Number of probes: 60x6=360
> Probe name prefix: Cg_consensus
> 
> Probe set 5: Cg genome Pho4 nonconsensus binding sites
> Number of probes: 90x6=540
> Probe name prefix: Cg_nonconsensus
> 
> Probe set 6: BLI library (UASp2 context) with random flanks
> Number of probes: 100x6=600
> Probe name prefix: BLI_UASp2
> 
> Probe set 7: Cg genome negative controls (regions not bound)
> Number of probes:150x6=900
> Probe name prefix: Cg_neg
> 
> Thanks,
> Lindsey

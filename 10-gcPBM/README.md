_Questions_

- How can CgPho4 bind to 50% more genomic sites with a E-box motif than ScPho4 does in the same genome?

_Hypotheses_

- CgPho4 is more capable of competing with Cbf1 (or other competitors) at the genomic sites where ScPho4 is normally outcompeted, implying that it has higher affinity for those sites with unpreferred flanking base pairs for ScPho4.
- The additional sites have higher nucleosome occupancy and CgPho4 is more capable of accessing nucleosome-occluded sites -- either by outcompeting the nucleosomes or by direct binding to those sites.
- The additional genomic sites don’t encode Pho2 binding motifs; CgPho4 doesn’t require the help of Pho2 and thus only it can bind to those sites.

_Approach_

What can gcPBM measure?

- Binding specificity differences between paralogs or TFs in the same family with highly similar core motifs (able to interrogate a wider range of base pairs).
Limitation of the gcPBM

- Cannot directly compare the log intensity between two chips (active [TF]) can be different -- developed a WLSR approach to use replicate chips to learn the variance and use that to determine the probes that are differentially bound by the paralogous TFs.


_Experiments_

- Use gcPBM to obtain the individual Pho4 ortholog’s binding energy landscape (extended motif in 36 bp context) and compare the preference between ScPho4 and CgPho4 (we are also interested in comparing more distantly related orthologs).
- Use “competition PBM” to determine if one ortholog has higher affinity than the other at a broad group of target sites.


_Reference_

    Shen, Ning, Jingkang Zhao, Joshua L. Schipper, Yuning Zhang, Tristan Bepler, Dan Leehr, John Bradley, John Horton, Hilmar Lapp, and Raluca Gordan. “Divergence in DNA Specificity among Paralogous Transcription Factors Contributes to Their Differential In Vivo Binding.” Cell Systems 6, no. 4 (April 25, 2018): 470-483.e8. https://doi.org/10.1016/j.cels.2018.02.009.

_Other notes [HB]_

- My previous ChIP-derived motifs seem to suggest highly similar sequence preferences between ScPho4 and CgPho4. I should have, however, tried to learn a motif based on the CgPho4-only sites.
- Is there any hint that the two orthologs may differ in their flanking site preference?
- If there is indeed a difference, what would that mean? Will it help explain the expanded binding sites for CgPho4?
- One hypothesis is that CgPho4 is able to bind the E-box motifs that are either blocked by nucleosomes or don’t have strong Pho2 motifs around. The former seems a real possibility and we can test that hypothesis by analyzing published nucleosome occupancy data or generating our own under phosphate starvation. The second possibility is less well supported by existing evidence: Pho2 seems to be the “follower”, going where the main TF goes and not too picky about its own motif. That could just be my wrong intuition though.
- If the two orthologs don’t differ in their preference, could they still differ in affinity? Would that be detectable through the gcPBM? Should be, I think!
- At least we should be able to measure the two Pho4’s binding to those exact sites without the nucleosome or help from Pho2. This should help us rule out some possibilities?
- In Yuning’s experiment, the array had 1500 probes, which include almost all CACGTG instances along with a non-specific probe set and a MITOMI-based Kd calibration set. Because they were very lenient when choosing the probes based on the ChIP data, the existing probe set almost certainly includes the ChIP-peaks in my 2017 paper for both ScPho4 and CgPho4.
- Raluca and Yuning both suggested adding the C. glabrata genome probes -- in that experiment I found a lot more binding peaks without a consensus motif.
- The probes are 36 bp long, with an additional 24 bp that encode the complementary DNA
- When doing protein purification, they were not too concerned about removing the degradation product, as long as they don’t contain the DBD. For ScPho4, they have been using frozen aliquots even 1-2 years after the initial purification. Note that compared with BLi, which is very sensitive in revealing protein quality issues, techniques like EMSA are a lot more tolerant because what’s being measured is the fraction of DNA that is shifted.
- PBM doesn’t give absolute binding affinities. However, they did include a Kd calibration probe set on the array, and they found that as long as the protein is not saturating the array, the log intensity measurements scale linearly with log Kd.
- They have used PBM to measure cooperative binding for Ets-1 and Runt. We should think about this for Pho4 and Pho2.

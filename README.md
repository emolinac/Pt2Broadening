# Pt2 Broadening
Various codes related to the Pt2 Broadening analysis.

## Organization
- essentials : Stores all the codes related to the operations previous to visualization of the broadening
- utilities : Stores codes not related to broadening or root files that are used by *essentials* and *visualization*
- visualization : Stores all the codes that are used to visualize the broadening in function of different kinematical variables

## Quick Tutorial
To obtain the broadening we need:
- A root file that contains the PhiPQ distributions in Q2, Nu, Zh and Pt2 bins
- A root file that contains an Ntuple with the binning (*check utilities*)

The steps are:
1. Integrate the PhiPQ distributions with *MACRO_Integration_Phi.cpp*
2. Treat the Pt2 distributions with *MACRO_Pt2_processing.cpp*

After the previoous steps, there are multiple choices about what to do. There are MACROS for every choice!
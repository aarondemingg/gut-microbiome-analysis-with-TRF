# gut-microbiome-analysis-with-TRF
 A project which I worked on my PhD. some of the codes which I used to process the data. The paper is currently undergoing a review for publication


# Introduction to the analysis I did
## A small background on my topic

Tocotrienol-Rich Fraction (TRF) is a form of vitamin E that is derived from sources like palm oil, rice bran, and annatto seeds. Unlike the more common tocopherols, tocotrienols have unique chemical structures that provide them with distinct biological properties. TRF has been shown to possess potent antioxidant and anti-inflammatory properties. These properties can help modulate immune responses, potentially enhancing the body's ability to fight off infections and reducing inflammation. TRF can influence the gut microbiome by promoting a healthier balance of gut bacteria. This, in turn, can positively affect overall health and immune function. The antioxidant properties of TRF can also help protect the gut lining from oxidative stress and inflammation. In summary, TRF has the potential to modulate immune responses and improve gut health through its antioxidant and anti-inflammatory properties. 

For this project, I investigated the role of TRF on the gut microbiome. As this is the first time I am doing this, please forgive me for any mistakes or inconsistencies that will be made in this project.

# availabilty of data

I will make files available once the study has been published

## Installation

You will need to install [R](https://www.r-project.org/) and an IDE of your choice (Mine's [Rstudio](https://posit.co/download/rstudio-desktop/)) for this project.

# workflow

I first cleaned and processed the raw data with [DADA2](https://benjjneb.github.io/dada2/tutorial.html) followed by [Paprica](https://github.com/bowmanjeffs/paprica) see here. The full data analysis [here](https://github.com/aarondemingg/gut-microbiome-analysis-with-TRF/blob/main/R/analysis.r)

# References 

Bowman, Jeff S., and Hugh W. Ducklow, 2015. Microbial Communities Can Be Described by Metabolic Structure: A General Framework and Application to a Seasonally Variable, Depth-Stratified Microbial Community from the Coastal West Antarctic Peninsula. PloS one, 10.8, e0135868

Callahan, B.J., McMurdie, P.J., Rosen, M.J., Han, A.W., Johnson, A.J.A. and Holmes, S.P., 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), pp.581-583.

Lahti, L., Sudarshan, S. and Ernst, F.M., 2021. Orchestrating Microbiome Analysis with Bioconductor [Beta Version](https://microbiome.github.io/OMA/docs/devel/#ref-OMA) 

McMurdie, P.J. and Holmes, S., 2013. phyloseq: an R package for reproducible interactive analysis and graphics of microbiome census data. PloS one, 8(4), p.e61217.
# wodu-bbwd
Code and data to use eBird S&amp;T products to analyze the effect of black-bellied whistling duck expansion on wood duck population trend

Stillman et al. 2025 "A Framework for Assessing the Habitat Correlates of Spatially Explicit Population Trends" used interpretable machine learning generalized additive model (GAM) to understand the relationships between land cover and spatially explicit bird population trends. Their analysis described a method to use on land cover variables from eBird Status and Trends products and how those effect species trends. 

This project aims to extend that modeling framework to include abundance and trend of another species as a variable in the GAM analysis. Abundance and trend data from other species in eBird S&T products are available at the same scale using the same statistical methods.

Black-bellied whistling ducks have expanded their range in the United States over the past 10-20 years, particurlarly in the southeastern US. Their range has expanded into areas where breeding wood ducks are also present, and because these are both cavity nesting species there is potentail for competition. There are no long-term surveys of either species in the southeastern US outside of the Breeding Bird Survey, a survey not well-suited for wetland-dependant species. This project can leverage the broad spatial extend of eBird S&T products to simultaniously estimate the effect of landscape variables and the changing distribution of other species on wood duck population trends.

In addition to habitat variables used in eBird S&T products, there may be more localized data sources to better estimate the effect of changing landcover on species trend. Forest Inventory Analysis data is available for the US, and may offer useful data pertinent to these species.

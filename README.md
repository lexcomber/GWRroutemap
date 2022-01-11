# A route map for successful applications of Geographically Weighted Regression

Alexis Comber<sup>1*</sup>, Chris Brunsdon<sup>2</sup>, Martin Charlton<sup>2</sup>, Guanpeng Dong<sup>3</sup>, Rich Harris<sup>4</sup>, Binbin Lu<sup>5</sup>, Yihe Lü<sup>6</sup>, Daisuke Murakami<sup>7</sup>, Tomoki Nakaya<sup>8</sup>, Yunqiang Wang<sup>9</sup>, Paul Harris<sup>10</sup>

<sup>1</sup> School of Geography, University of Leeds, Leeds, UK.\
<sup>2</sup> National Centre for Geocomputation, Maynooth University, Maynooth, Ireland.\
<sup>3</sup> Key Research Institute of Yellow River Civilization and Sustainable Development, Henan University, Kaifeng, China.\
<sup>4</sup> School of Geographical Sciences, University of Bristol, Bristol, UK.\
<sup>5</sup> School of Remote Sensing and Information Engineering, Wuhan University, Wuhan, China.\
<sup>6</sup> State Key Laboratory of Urban and Regional Ecology, Research Center for Eco-Environmental Sciences, Chinese Academy of Sciences; Joint Center for Global Change Studies; University of Chinese Academy of Sciences, Beijing, China.\
<sup>7</sup> Department of Statistical Data Science, Institute of Statistical Mathematics, Tachikawa, Japan.\
<sup>8</sup> Graduate School of Environmental Studies, Tohoku University, Sendai, Japan.\
<sup>9</sup> State Key Laboratory of Loess and Quaternary Geology, Institute of Earth Environment, Chinese Academy of Sciences, Xi’an, China.\
<sup>10</sup> Sustainable Agriculture Sciences North Wyke, Rothamsted Research, Okehampton, UK.

<sup>*</sup> contact author: a.comber@leeds.ac.uk

## Abstract

Geographically Weighted Regression (GWR) is increasingly used in spatial analyses of social and environmental data. It allows spatial heterogeneities in processes and relationships to be investigated through a series of local regression models rather than a single global one. Standard GWR assumes that relationships between the response and predictor variables operate at the same spatial scale, which is frequently not the case. To address this, several GWR variants have been proposed. This paper describes a route map to decide whether to use a GWR model or not, and if so which of three core variants to apply: a standard GWR, a mixed GWR or a multiscale GWR (MS-GWR). The route map comprises 3 primary steps that should always be undertaken: (1) a basic linear regression, (2) a MS-GWR, and (3) investigations of the results of these in order to decide whether to use a GWR approach, and if so for determining the appropriate GWR variant. The paper also highlights the importance of investigating a number of secondary issues at global and local scales including collinearity, the influence of outliers, and dependent error terms. Code and data for the case study used to illustrate the route map are provided.

**Keywords**: Spatially varying coefficient model; non-stationarity; spatial heterogeneity; autocorrelation; regression

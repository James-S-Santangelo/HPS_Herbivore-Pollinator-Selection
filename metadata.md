# Metadata for datasets
##### Manuscript: Herbivores and plant defenses affect selection on plant reproductive traits more strongly than pollinators
##### Authors: James S. Santangelo, Ken A. Thompson and Marc T.J. Johnson
##### Journal: Journal of Evolutionary Biology

### Data measured on all plants. 'datExp' in R script.

| Column | Description | Type |
|--------|-------------|------|
| Plant  | ID for plant. Use Genotype instead. | Integer |
| Parent.Genotype | Unique plant genotype | Factor |
| Stolon | Genotype replicate. Range: 1-16 | Integer |
| Label | Label of plant at field site | Factor |
| Row | Row position of plant at field site. Range: 1-23 | Integer |
| Column | Column position of plant at field site. Range: 1-36 | Integer |
| Pollination | Level of pollination treatment. "Supp" or "Open" | Factor |
| Herbivory | Level of herbivory treatment. "Ambient" or "Reduced" | Factor |
| HCN | Presence (1) or absence (0) of HCN. | Integer |
| Glycosides.Ac | Presence (1) or absence (0) of cyanogenic glucosides. | Integer |
| Linamarase.Li| Presence (1) or absence (0) of Linamarase. | Integer |
| Prop.mass | Mass (g) of stolon planted before experiment | Float |
| Flower.date | Number of days to first flower | Integer |
| Mammal.herb | Presence (1) or absence (0) of damage by vertebrate herbivores | Integer |
| Status | "ALIVE" or "DEAD" | Factor |
| Leaf.avg.dmg.1 | Mean % leaf area lost during first herbivore damage survey | Float |
| Leaf.avg.dmg.2 | Mean % leaf area lost during second herbivore damage survey | Float |
| Leaf.avg.dmg.3 | Mean % leaf area lost during third herbivore damage survey | Float |
| Avg.Bnr.Wdth| Average width (mm) of banner petals | Float |
| Avg.Bnr.Ht | Average length (mm) of banner petals | Float |
| Avg.bnr.dmg | Average % banner petal lost due to herbivore damage | Float |
| Bag.Inf.Wht | Mass (g) bagged inflorescenes | Float |
| Bag.Seed.Num | Count of seeds from bagged inflorescences | Integer |
| Bag.Seed.Wht | Madd (g) of seeds fromg bagged inflorescneces | Float |
| Biomass.plant | Biomass (g) of vegetative plant tissue | Float |
| Num.flwrs | Mean number of flowers per inflorescence | Integer |
| Total.Seed.mass| Total mass (g) of seeds across all inflorescences | Float |
| Total.Inf  | Total number of inflorescences produced by plant | Integer |
| Seeds.Inf | Mean mass (g) number of seeds per inflorescence | Float |
| Seed.only | Mass (g) of only seeds. Use Total.Seed.mass instead | Float |
| Block | Spatial Block. Range: A-F. | Factor |

### Data for trial assessing the effects of pesticides on plant fitness. 'datIns' in R script.

| Column | Description | Type |
|--------|-------------|------|
| Date.planted  | Date stolon was planted | Date string |
| Parent.plant | Unique plant ID | Integer |
| Stolon | Genotype replicate. Range: 1-30 | Integer |
| Label | Label of plant | Integer |
| Tray | Tray in which plant was placed. Range: 1-12 | Integer |
| Position | Position of plant in tray. Range: 1-8 | Integer |
| Treatment | Concatenation of 'Molluscicide' and 'Insecticide' treatment columns | Factor |
| Molluscicide | Level of molluscicide treatment. "No" (absent) or "Yes" (present) | Factor |
| Insecticide | Level of insecticide treatment. "None", "Low", or "High" | Factor |
| Prop.mass | Mass (g) of stolon planted before experiment | Float |
| In.Leaves | Initial number of leaves | Integer |
| Num.Flowers | Number of flowers of pollinated inflorescence | Integer |
| Num.Seeds | Number of seeds produced by pollinated inflorescence | Integer |
| Seeds.flw | Number of seeds per flower | Float |
| Biomass | Final mass (g) of vegetative plant tissue | Float |

### Data for pollinator observations during experiment. 'datPoll.obs' in R script.

| Column | Description | Type |
|--------|-------------|------|
| Date  | Date on which observations were carried out | Date string |
| Pollinator | Sequential number given to each pollinator theat was tracked | Integer |
| Group | Pollinator morphogroup | Factor |
| Plant | Label of plant | Factor |
| Insecticide | Level of pesticide treatment. "Yes" (treated) or "No" (untreated) | Factor |
| Num.Inf | Number of inflorescences on plant at time of observation. | Integer |
| Time.s | Time (s) spent foraging by pollinator | Float |

### Data for genotype means generated from the data measured on all plants. Used for genotypic selection analysis to examine effects of treatments on patterns of selection. 'GTSelnData' in R script.

| Column | Description | Type |
| ------ | ----------- | ---- |
| Genotype | Plant genotype | Factor |
| HCN | Presence (1) or absence (0) of HCN. | Integer |
| Herbivory| Level of herbivory treatment. "Ambient" or "Reduced" | Factor |
| Pollination| Level of pollination treatment. "Supp" or "Open" | Factor |
| Glycosides.Ac| Presence (1) or absence (0) of cyanogenic glucosides. | Integer |
| Linamarase.Li| Presence (1) or absence (0) of Linamarase. | Integer |
| Flwr.date | Mean number of days to first flower | Integer |
| Bnr.wdth | Mean width (mm) of banner petals | Float |
| Bnr.ht | Mean length (mm) of banner petals | Float |
| Biomass | Mean biomass (g) of vegetative plant tissue | Float |
| Infl | Mean number of inflorescences produced by plant | Integer |
| Flwrs | Mean number of flowers per inflorescence | Integer |
| AF.Seeds | Mean total mass (g) of seeds | Float |
| RF.Seed | Relative fitness | Float |
| Flwr.date.T | Transformed date to first flower | Float
| Bnr.wdth.T | Transformed banner petal width | Float |
| Bnr.ht.T | Transformed banner petal length | Float |
| Biomass.T | Transformed plant biomass | Float |
| Infl.T | Tranformed number of inflorescences | Float |
| Flwrs.T | Transformed number of flowers per inflorescnece | Float |
| Flwr.date.S | Standardized date to first flower | Float |
| Bnr.wdth.S | Standardized banner petal width | Float |
| Bnr.ht.S | Standardized  banner petal length | Float |
| Biomass.S | Standardized plant biomass | Float |
| Infl.S | Standardized number of inflorescences | Float |
| Flwrs.S | Standardized number of flowers per inflorescences | Float |
| AF.Seeds.S | Standardized absolute fitness. Not used. | Float |

### Data for genotype means generated from the data measured on all plants. Used for genotypic selection analysis to examine effects of vole damage on patterns of selection. 'GTSelnData' in R script.

| Column | Description | Type |
| ------ | ----------- | ---- |
| Genotype | Plant genotype | Factor |
| Mammal.herb | Presence (1) or absence (0) of damage by vertebrate herbivores | Integer |
| HCN | Presence (1) or absence (0) of HCN. | Integer |
| Flwr.date | Mean number of days to first flower | Integer |
| Bnr.wdth | Mean width (mm) of banner petals | Float |
| Bnr.ht | Mean length (mm) of banner petals | Float |
| Biomass | Mean biomass (g) of vegetative plant tissue | Float |
| Infl | Mean number of inflorescences produced by plant | Integer |
| Flwrs | Mean number of flowers per inflorescence | Integer |
| AF.Seeds | Mean total mass (g) of seeds | Float |
| RF.Seed | Relative fitness | Float |
| Flwr.date.T | Transformed date to first flower | Float
| Bnr.wdth.T | Transformed banner petal width | Float |
| Bnr.ht.T | Transformed banner petal length | Float |
| Biomass.T | Transformed plant biomass | Float |
| Infl.T | Tranformed number of inflorescences | Float |
| Flwrs.T | Transformed number of flowers per inflorescnece | Float |
| Flwr.date.S | Standardized date to first flower | Float |
| Bnr.wdth.S | Standardized banner petal width | Float |
| Bnr.ht.S | Standardized  banner petal length | Float |
| Biomass.S | Standardized plant biomass | Float |
| Infl.S | Standardized number of inflorescences | Float |
| Flwrs.S | Standardized number of flowers per inflorescences | Float |






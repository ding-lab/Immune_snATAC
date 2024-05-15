# This folder contains script that is used to calculate AUCell scores for all cells (not only those that were used for regulons inference).

1. ```Make_loom_obj_all_cells.20231127.R``` -- First, make a loom object that contains all cells.

2. ```AUCell_new_regulons_All_cells.20231127.R``` -- Then, use all cells from the object to quantify AUCell scores for the regulons identified after filtering by their frequency.
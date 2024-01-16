For estimating the retention/loss patterns of the rho-derived duplicates in different subfamilies, custom scripts were designed to detect gene trees in large-scale.

These scripts require running in Mac with the arm64 CPU architecture (the Python version=3.9.6).

Below is a step-by-step guide for analyzing a small test dataset.

Step 1, gene trees will be refined by the removal of long branches using the following command.
python3 TreeAnalyses_01_Prune_long_branches.py TestTree.list

Step2, gene trees will be rooted using the following command.
python3 TreeAnalyses_02_root_using_EvoMADs.py HOG.pineapple ref.topology.tree Grass_tips
[In this step, our program invoked the minimal ancestor deviation approach (Tria et al., 2017 Nature Ecology & Evolution; see their mad.py script).]

Step3, gene trees will be reconciled to species-tree for estimating retention/loss patterns.
python3 TreeAnalyses_03_estimate_retention.py Rho.tree_list ref.topology.tree Grass_tips Rho.tree_list.dup.txt

Step4, filter by gene pairs mapped at the nodes earlier than Poaceae.
python3 TreeAnalyses_04_annotation.py
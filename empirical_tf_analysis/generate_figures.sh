python3 data_generation.py forbes/scores.txt "Galaxy52-[GSuite_(29_-_Filter_categorical_gSuite)_-_ready_for_analysis].gsuite" data/forbes.csv
python3 data_generation.py jaccard/scores.txt "Galaxy52-[GSuite_(29_-_Filter_categorical_gSuite)_-_ready_for_analysis].gsuite" data/jaccard.csv 
python3 data_generation.py tetra/scores.txt "Galaxy52-[GSuite_(29_-_Filter_categorical_gSuite)_-_ready_for_analysis].gsuite" data/tetra.csv 

Rscript analyze_data_v2.r 

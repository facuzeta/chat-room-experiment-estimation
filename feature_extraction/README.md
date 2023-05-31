## About:
This folder contains the code to extract features from the chat logs of the experiment.

## Files:

### compute_metrics_and_create_csvs.py:

This script takes one raw data file and generates the output files with the computed metrics.

For example:

```bash
python3 compute_metrics_and_create_csvs.py ../raw_data/data_experiment_WITHOUT_BAIAS.json
python3 compute_metrics_and_create_csvs.py ../raw_data/data_experiment_WITH_BAIAS_FERMI.json
python3 compute_metrics_and_create_csvs.py ../raw_data/data_experiment_WITH_BAIAS_NUMBERS.json
python3 compute_metrics_and_create_csvs.py ../raw_data/data_experiment_WITH_BAIS_FERMI_INDIVIDUAL.json
```

Generates:

```bash
data_experiment_WITHOUT_BAIAS__by_participant.csv
data_experiment_WITHOUT_BAIAS__with_features.csv
data_experiment_WITH_BAIAS_FERMI__by_participant.csv
data_experiment_WITH_BAIAS_FERMI__with_features.csv
data_experiment_WITH_BAIAS_NUMBERS__by_participant.csv
data_experiment_WITH_BAIAS_NUMBERS__with_features.csv
data_experiment_WITH_BAIS_FERMI_INDIVIDUAL__by_participant.csv
data_experiment_WITH_BAIS_FERMI_INDIVIDUAL__with_features.csv
```

### metrics.py:

This file contains the functions to compute the metrics.

### dic_word_vec.json:
This file contains the word embeddings for the words in the chat logs using FastText. This file is generated by the script `create_dic_word_vec.py`. The purpose of this file is to avoid to load Fasttext every time we want to compute the metrics becasue this model requieres a lot of memory.

### create_dic_word_vec.py:
This file generates the file `dic_word_vec.json` using FastText.
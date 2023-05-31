import json
from metrics import (
    compute_fermi_whitelist,
    compute_exchange_numbers,
    compute_ff,
)
import pandas as pd
import numpy as np
import json
from tqdm import tqdm
import sys


EXPERIMENT_WITHOUT_INTERVENTION = 1
EXPERIMENT_INTERVENTION_NUMBERS = 3
EXPERIMENT_INTERVENTION_FERMI = 4
EXPERIMENT_INTERVENTION_FERMI_INDIVIDUAL = 5

DIC_EXP_LABELS = {}
DIC_EXP_LABELS[EXPERIMENT_WITHOUT_INTERVENTION] = "Exp. without intervention"
DIC_EXP_LABELS[EXPERIMENT_INTERVENTION_NUMBERS] = "Exp. with intervention (numbers)"
DIC_EXP_LABELS[EXPERIMENT_INTERVENTION_FERMI] = "Exp. with intervention (Fermi)"
DIC_EXP_LABELS[EXPERIMENT_INTERVENTION_FERMI_INDIVIDUAL] = "Exp. with intervention (Fermi) individual"


def add_features_to_data(data):
    """
    Add features to data: for every group and every chat, compute scores
    """

    for d in tqdm(data):
        for stage_id in d["stage_2"]:
            features = {}
            features.update(compute_fermi_whitelist(d["stage_2"][stage_id]))
            features.update(compute_exchange_numbers(
                d["stage_2"][stage_id], d["stage_1"]))
            features.update(compute_ff(d["stage_2"][stage_id], dic_word_vec))
            d["stage_2"][stage_id]["features"] = features

    res = []
    for d in tqdm(data):
        for stage_2_i in d["stage_2"]:
            r = {}
            r["group_id"] = d["group_id"]
            r["experiment_id"] = d["experiment_id"]
            r["experiment_label"] = DIC_EXP_LABELS[r["experiment_id"]]

            r["question_id"] = d["stage_2"][stage_2_i]["question_id"]

            r["raters_fermi_mean"] = np.mean(
                d["stage_2"][stage_2_i]["raters"]["FERMI_VALUATION"])
            r["raters_number_mean"] = np.mean(
                d["stage_2"][stage_2_i]["raters"]["NUMBERS_VALUATION"])

            # mean answers all subject in s1 for the same question
            r["stage_2_i"] = stage_2_i
            r["s1_mean_value"] = np.nanmean(
                [e["answer_value"] for e in d["stage_1"] if e["question_id"]
                    == d["stage_2"][stage_2_i]["question_id"]]
            )
            r["consensus_value"] = d["stage_2"][stage_2_i]["consensus_value"]
            r.update(dict(d["stage_2"][stage_2_i]["features"]))
            res.append(r)
    df = pd.DataFrame(res)

    # remove data with no consensus
    df = df[~(df.consensus_value.isna())]
    return df


def compute_metric_by_participant(data):
    """
    Compute metrics by participant
    """

    def if_is_not_defined_return_0(d, list_keys):
        try:
            dd = d
            for k in list_keys:
                dd = dd[k]
            return dd
        except:
            return 0
    new_res = []
    for d in tqdm(data):

        participants = sorted(set([e["participant_id"] for e in d["stage_1"]]))
        for p in participants:
            for question_id in sorted(set([e['question_id'] for e in d['stage_1']])):
                new_r = dict()
                new_r["subject_id"] = p
                new_r["group_id"] = d["group_id"]

                new_r["question_id"] = question_id
                new_r["experiment_id"] = d["experiment_id"]
                new_r["experiment_label"] = DIC_EXP_LABELS[d["experiment_id"]]


                possible_ans_s1 = [
                    r for r in d["stage_1"] if r["participant_id"] == p and r["question_id"] == new_r["question_id"]
                ]
                if len(possible_ans_s1) > 0:
                    new_r["s1_value"] = possible_ans_s1[0]["answer_value"]
                    new_r["s1_confidence"] = possible_ans_s1[0]["confidence_value"]

                possible_ans_s3 = [
                    r for r in d["stage_3"] if r["participant_id"] == p and r["question_id"] == new_r["question_id"]
                ]
                if len(possible_ans_s3) > 0:
                    new_r["s3_value"] = possible_ans_s3[0]["answer_value"]
                    new_r["s3_confidence"] = possible_ans_s3[0]["confidence_value"]

                # check if participant discussed this question
                dic_question_id_stage2_i = {
                    dd['question_id']: s2_i for s2_i, dd in d['stage_2'].items()}
                if question_id in dic_question_id_stage2_i:
                    s2_i = dic_question_id_stage2_i[question_id]
                    # add some features
                    new_r["n_fermi_whitelist"] = if_is_not_defined_return_0(
                        d, ["stage_2",s2_i,"features","features_fermi_whitelist_count",p]
                    )
                    new_r["n_number_en_s1"] = if_is_not_defined_return_0(
                        d, ["stage_2",s2_i,"features","features_number_whitelist_count",p]
                    )
                    new_r["n_interventions"] = len(
                        [e for e in d["stage_2"][s2_i]["chat"] if e["participant_id"] == p])
                    new_r["n_words"] = sum(
                        [1 + e["text"].count(" ") for e in d["stage_2"]
                         [s2_i]["chat"] if e["participant_id"] == p]
                    )

                new_res.append(dict(new_r))
    return pd.DataFrame(new_res)


if __name__ == "__main__":
    # load word vectors
    dic_word_vec = json.load(open('dic_word_vec.json'))

    # use first argument as input file
    fn = sys.argv[1]
    print("Using file", fn)

    # read data
    with open(fn, "r") as f:
        data = json.load(f)
    print("Read", len(data), "groups")

    # add features to data
    df = add_features_to_data(data)

    # save data groups
    fnout = fn.replace(".json", "__with_features.csv").replace(
        "raw_data", "output").split("/")[-1]
    print("Saving data in", fnout)
    features_to_export = [  'group_id',
                            'experiment_id',
                            'experiment_label',
                            'question_id',
                            'stage_2_i',
                            's1_mean_value',
                            'consensus_value',
                            'raters_fermi_mean',
                            'raters_number_mean',
                            'features_fermi_whitelist',
                            'features_number_whitelist',
                            'features_ff_mean',
                            'features_ff_std',
                            'features_ff_min',
                            'features_ff_25%',
                            'features_ff_50%',
                            'features_ff_75%',
                            'features_ff_max',
                            'features_ff_percentile_90']

    # remove chats without consensus
    df = df[df.consensus_value != -1]

    df[features_to_export].to_csv(fnout, index=False)

    # save data by participant
    df_participant = compute_metric_by_participant(data)
    features_to_export = [
        "subject_id",
        "group_id",
        "experiment_id",
        "experiment_label",
        "question_id",
        "s1_value",
        "s1_confidence",
        "s3_value",
        "s3_confidence",
        "n_number_en_s1",
        "n_fermi_whitelist",
        "n_interventions",
        "n_words",
    ]
    fnout = fn.replace(".json", "__by_participant.csv").replace(
        "raw_data", "output").split("/")[-1]
    df_participant[features_to_export].to_csv(fnout, index=False)

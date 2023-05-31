import json
import fasttext
import numpy as np
import os


# replace this with the correct filename for the fasttext model
# The model is available in https://fasttext.cc/docs/en/crawl-vectors.html
FASTTEX_SPANISH_MODEL_FN = './cc.es.300.bin'
PATH = '../raw_data'

def is_only_letter(txt):
    "Check if a string is composed only by letters"
    return len(txt) > 0 and all([c.isalpha() for c in txt])


def get_norm_vector(ft, w):
    "Get the normalized vector for a word"
    v = ft.get_word_vector(w)
    return (v/np.linalg.norm(v)).tolist()


def create_dic_word_vec():
    "Create a dictionary with the vector for each word in the dataset"

    # read raw data
    data = []
    fns = os.listdir(PATH)
    for fn in fns:
        fin = open(os.path.join(PATH, fn), 'r')
        data += json.load(fin)
        fin.close()

    # list all words in dataset
    # participant answers
    words = sorted(set(([w.lower() for g in data for ll in g['stage_2'].values(
    ) for r in ll['chat'] for w in r['text'].split() if is_only_letter(w)])))
    # questions
    words += [w.lower() for w in '¿Cuántos países tienen territorio en el hemisferio sur, excluyendo el continente antártico? ¿Cuántos capítulos tiene la serie estadounidense de televisión "Friends"? ¿Cuántos Kilos de helado come por año una persona en Argentina? Hasta 2021, ¿en cuántas finales de la Copa Libertadores participó un equipo argentino? ¿Cuál es la longitud en metros del recorrido total de la línea A? (subte de CABA) ¿Cuántas personas (en MILLONES) viven en América del Sur? ¿Cuántos escalones tiene la escalera interna del obelisco? ¿Cuántos goles se anotaron en total en el mundial Brasil 2014?'.split()]
    words = sorted(set(words))

    # load fasttet model
    ft = fasttext.load_model(FASTTEX_SPANISH_MODEL_FN)

    # save dic_word_vec.json
    dic_word_vector = {w: get_norm_vector(ft, w) for w in words}
    fout = open('dic_word_vec.json', 'w')
    json.dump(dic_word_vector, fout)
    fout.close()

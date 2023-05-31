import json
import re

import numpy as np
import pandas as pd
from tqdm import tqdm
from nltk.corpus import stopwords

FEATURES_FERMI_WHITELIST_NAME = "features_fermi_whitelist"
FEATURES_NUMBER = "features_number_whitelist"
FEATURES_FF = "features_ff"


def compute_fermi_whitelist(stage_2_i):
    # compute fermi score using whitelist method
    FERMI_WHITELIST_DICT = {
        1: [
            "partido",
            "partidos",
            "mundial",
            "grupo",
            "grupos",
            "64",
            "pais",
            "paíspaises",
            "final",
            "cuartos",
            "octavos",
            "semi",
            "semifinal",
            "semifinales",
            "promedio",
            "penales",
            "penal",
        ],
        2: ["metro", "metros", "piso", "pisos", "edificio", "edificios", "escalon"],
        3: [
            "pais",
            "país",
            "argentina",
            "bolivia",
            "brasil",
            "brasíl",
            "chile",
            "colombia",
            "ecuador",
            "guayana",
            "paraguay",
            "peru",
            "perú",
            "surinam",
            "uruguay",
            "venezuela",
        ],
        4: [
            "cuadra",
            "cuadras",
            "metro",
            "metros",
            "estacion",
            "estaciones",
            "calle",
            "calles",
            "km",
            "entre",
            "larga",
            "plaza de mayo",
            "peru",
            "perú",
            "piedras",
            "lima",
            "saenz peña",
            "congreso",
            "pasco",
            "alberti",
            "plaza miserere",
            "loria",
            "castro barros",
            "rio de janeiro",
            "acoyte",
            "primera junta",
        ],
        5: [
            "river",
            "boca",
            "racing",
            "independiente",
            "futbol",
            "campeonato",
            "compeonatos",
            "velez",
            "ganadas",
            "gano",
            "san lorenzo",
        ],
        6: ["dia", "semana", "mes", "año", "1/4", "come", "kilo", "cuarto", "medio", "1/2", "cuartos"],
        7: ["temporadas", "temporada", "capitulos", "capitulo", "años", "año"],
        8: [
            "america",
            "américa",
            "áfrica",
            "africa",
            "oceania",
            "oceanía",
            "norte",
            "sur",
            "Angola",
            "Argentina",
            "Australia",
            "Bolivia",
            "Botsuana",
            "Brasil",
            "Burundi",
            "Chile",
            "Colombia",
            "Comoras",
            "Ecuador",
            "Fiyi",
            "Gabón",
            "Guinea",
            "Ecuatorial",
            "Indonesia",
            "Islas",
            "Salomón",
            "Kenia",
            "Kiribati",
            "Lesoto",
            "Madagascar",
            "Malaui",
            "Maldivas",
            "Mauricio",
            "Mozambique",
            "Namibia",
            "Nauru",
            "Nueva Zelanda",
            "Papúa",
            "Nueva",
            "Guinea",
            "Paraguay",
            "Perú",
            "República Congo",
            "República Democrática del Congo",
            "Ruanda",
            "Samoa",
            "Santo Tomé y Príncipe Seychelles",
            "Somalia",
            "Suazilandia",
            "Sudáfrica",
            "Tanzania",
            "Timor Oriental",
            "Tonga",
            "Tuvalu",
            "Uganda",
            "Uruguay",
            "Vanuatu",
            "Zambia",
            "Zimbabue",
        ],
    }

    chat_r = stage_2_i["chat"]
    question_id = stage_2_i["question_id"]
    chat = pd.DataFrame(chat_r)
    if len(chat) == 0:
        return {}

    r = {}
    r[FEATURES_FERMI_WHITELIST_NAME] = np.mean(
        [any((c.lower() in t.lower() for c in FERMI_WHITELIST_DICT[question_id]))
         for t in chat.text]
    )
    r[FEATURES_FERMI_WHITELIST_NAME + "_count"] = (
        pd.DataFrame(
            [
                {
                    "participant_id": rr.participant_id,
                    "count": sum([rr.text.lower().split(" ").count(c.lower())
                                  for c in FERMI_WHITELIST_DICT[question_id]]),
                }
                for rr in chat.iloc
            ]
        )
        .groupby("participant_id")
        .sum()
        .to_dict()["count"]
    )
    return r


# compute exchange nunmber score method
def compute_exchange_numbers(stage_2_i, stage_1):
    def only_words_or_letters(w):
        return re.sub(r"[^\w]+", "", w)

    chat_r = stage_2_i["chat"]
    question_id = stage_2_i["question_id"]
    chat = pd.DataFrame(chat_r)
    s1_answers = set([str(e["answer_value"])
                     for e in stage_1 if e["question_id"] == question_id])
    if len(chat) == 0:
        return {}
    r = {}
    r[FEATURES_NUMBER] = np.mean(
        [any((only_words_or_letters(w) in s1_answers for w in t.split(" ")))
         for t in chat.text]
    )
    r[FEATURES_NUMBER + "_count"] = (
        pd.DataFrame(
            [
                {
                    "participant_id": rr.participant_id,
                    "count": sum([only_words_or_letters(w) in s1_answers for w in rr.text.split(" ")]),
                }
                for rr in chat.iloc
            ]
        )
        .groupby("participant_id")
        .sum()
        .to_dict()["count"]
    )
    return r


def compute_ff(stage_2_i, dic_word_vec):
    "This function computes exchange number score method"
    def only_words_or_letters(w):
        return re.sub(r"[^\w ]+", "", w.lower())
    dic_questions_text = {'NUMBER_COUNTRIES_SOUTH_HEMISPHERE': '¿Cuántos países tienen territorio en el hemisferio sur, excluyendo el continente antártico?',
                          'NUMBER_CHAPTERS_FRIENDS_TV_SHOW': '¿Cuántos capítulos tiene la serie estadounidense de televisión "Friends"?',
                          'KILOGRAM_OF_ICECREAM_PER_PERSON_ARGENTINA': '¿Cuántos Kilos de helado come por año una persona en Argentina?',
                          'FINALS_LIBERTADORES_CUP_ARGENTINE_TEAMS': 'Hasta 2021, ¿en cuántas finales de la Copa Libertadores participó un equipo argentino?',
                          'LENGTH_UNDERGROUND_LINE_A': '¿Cuál es la longitud en metros del recorrido total de la línea A? (subte de CABA)',
                          'POPULATION_SOUTH_AMERICA': '¿Cuántas personas (en MILLONES) viven en América del Sur?',
                          'STEPS_OBELISCO': '¿Cuántos escalones tiene la escalera interna del obelisco?',
                          'GOALS_WORDLCUP_2014': '¿Cuántos goles se anotaron en total en el mundial Brasil 2014?',
                          }

    STOPWORDS_ES = set(stopwords.words('spanish'))
    STOPWORDS_ES.update(set(['entonces', 'vamos', 'bueno', 'bien', 'claro', 'confianza', 'igual',
                             'hola', 'acuerdo', 'solo', 'seguro', 'verdad', 'cuanto', 'poner',
                             'listo', 'capaz', 'ahora', 'medio', 'tipo', 'respuesta',  'cuenta',
                             'tener', 'haber']))

    chat = stage_2_i["chat"]
    if len(chat) == 0:
        return {}

    words = [w for line in chat for w in only_words_or_letters(line['text']).split() if w not in STOPWORDS_ES and len(w) > 2]
    M_chat = np.array([dic_word_vec[w] for w in words if w in dic_word_vec])
    if len(M_chat) == 0:
        return {}
    # print(words)
    M_trigger = np.array([dic_word_vec[w] for w in dic_questions_text[stage_2_i['question_code']].lower().split() if w in dic_word_vec])
    sim_todos_contra_todos = np.dot(M_trigger, M_chat.T).flatten()
    r = {FEATURES_FF+"_"+k: v for k,
         v in pd.Series(sim_todos_contra_todos).describe().to_dict().items()}
    r[FEATURES_FF+'_percentile_90'] = np.percentile(sim_todos_contra_todos, 90)
    return r

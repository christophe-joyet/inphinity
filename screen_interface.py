from tkinter import *
from tkinter.filedialog import askopenfilename

import network_chart_organism_script as ncos
import network_chart_couples_script as nccs
import matrix_similarity_script as mss
import get_features_from_organism_script as ftrs
import display_graphic


def openpath(e:Entry):
    path = askopenfilename(title="Ouvrir une image",filetypes=[('all files','.*')])
    e.set(path)

# fenetre principale
fenetre = Tk()
fenetre.geometry('700x500')
fenetre.title("Functions")

# ================================================================
# network_chart_script
# ================================================================

# Cadre 1 
Cadre_1 = LabelFrame(fenetre, text="Network Chart Script", height=80, width=550, relief=GROOVE, labelanchor='n')

# organism_id_1_network_chart_script
organism_id_1_txt = Label(Cadre_1, text = 'Organism id :', width=20)
organism_id_1_value = IntVar() 
organism_id_1 = Entry(Cadre_1, textvariable=organism_id_1_value, width=20)

# checkbutton_network_chart_script
is_phage_1_value = IntVar()
is_phage_1 = Checkbutton(Cadre_1, text="Phage", variable=is_phage_1_value)

# bouton network_chart_script
bouton_network_chart_script=Button(Cadre_1, text="Graphic : bacterium - phage", command= lambda: ncos.networkChartOrganismScript(organism_id=organism_id_1.get(), is_phage=is_phage_1_value.get()))

# display in the interface
Cadre_1.pack()
organism_id_1_txt.pack()
organism_id_1_txt.place(x=10, y=10)
organism_id_1.pack()
organism_id_1.place(x=150, y=10)
is_phage_1.pack()
is_phage_1.place(x=320, y=10)
bouton_network_chart_script.pack()
bouton_network_chart_script.place(x=10, y=35)

# ================================================================
# network_chart_couple_script
# ================================================================

# Cadre 2
Cadre_2 = LabelFrame(fenetre, text="Network Chart Couple Script", height=80, width=550, relief=GROOVE, labelanchor='n')

# liste lysis_couple
liste_lysis_couple_value = StringVar(Cadre_2)
liste_lysis_couple_value.set("Lysis Type")
liste_lysis_couple = OptionMenu(Cadre_2, liste_lysis_couple_value, "CLEAR_LYSIS", "SEMI_CLEAR_LYSIS", "OPAQUE_LYSIS", "CLEAR_LYSIS_1E7PLUS", 
                                                                    "SEMI_CLEAR_LYSIS_1E7PLUS", "CLEAR_LYSIS_1E7MINUS", "SEMI_CLEAR_LYSIS_1E7MINUS",
                                                                    "ALL_CLEAR_LYSIS", "ALL_SEMI_CLEAR_LYSIS")

# bouton network_chart_couple_script
bouton_network_chart_couple_script=Button(Cadre_2, text="graphic COUPLES: bacterium - phage", command= lambda: nccs.networkChartCouplesScript(parameter_type=liste_lysis_couple_value.get()))

# display in the interface
Cadre_2.pack()
bouton_network_chart_couple_script.pack()
bouton_network_chart_couple_script.place(x=10, y=10)
liste_lysis_couple.pack()
liste_lysis_couple.place(x=320, y=10)

# ================================================================
# create matrix similarity
# ================================================================

# Cadre 3
Cadre_3 = LabelFrame(fenetre, text="Create Matrix Similarity", height=125, width=550, relief=GROOVE, labelanchor='n')

# organism_id_2_network_chart_script
organism_id_2_txt = Label(Cadre_3, text = 'Organism id :', width=20)
organism_id_2_value = IntVar() 
organism_id_2 = Entry(Cadre_3, textvariable=organism_id_2_value, width=20)

# path
path_1_txt = Label(Cadre_3, text = 'Path to save :', width=20)
path_1_value = StringVar() 
path_1 = Entry(Cadre_3, textvariable=path_1_value, width=20)

# file_name_1
file_name_1_txt = Label(Cadre_3, text = 'Name of new file :', width=20)
file_name_1_value = StringVar() 
file_name_1 = Entry(Cadre_3, textvariable=file_name_1_value, width=20)

# bouton matrice
bouton_matrix=Button(Cadre_3, text="Create similitary matrix", command= lambda: mss.matrixSimilarityScript(file_name=file_name_1.get(), path=path_1.get(), organism_id=organism_id_2.get(), is_phage=is_phage_2_value.get()))

# checkbutton_network_chart_script
is_phage_2_value = IntVar()
is_phage_2 = Checkbutton(Cadre_3, text="Phage", variable=is_phage_2_value)

# display in the interface
Cadre_3.pack()
organism_id_2_txt.pack()
organism_id_2_txt.place(x=10, y=10)
organism_id_2.pack()
organism_id_2.place(x=150, y=10)
is_phage_2.pack()
is_phage_2.place(x=320, y=10)
path_1_txt.pack()
path_1_txt.place(x=10, y=30)
path_1.pack()
path_1.place(x=150, y=30)
file_name_1_txt.pack()
file_name_1_txt.place(x=10, y=50)
file_name_1.pack()
file_name_1.place(x=150, y=50)
bouton_matrix.pack()
bouton_matrix.place(x=10, y=75)

# ================================================================
# Extract features
# ================================================================

# Cadre 4
Cadre_4 = LabelFrame(fenetre, text="Extract features", height=120, width=550, relief=GROOVE, labelanchor='n')

# organism_id_1_network_chart_script
organism_id_3_txt = Label(Cadre_4, text = 'Bacterie id :', width=20)
organism_id_3_value = IntVar() 
organism_id_3 = Entry(Cadre_4, textvariable=organism_id_3_value, width=20)

# path
path_2_txt = Label(Cadre_4, text = 'Path to save :', width=20)
path_2_value = StringVar() 
path_2 = Entry(Cadre_4, textvariable=path_2_value, width=20)

# file_name_1
file_name_2_txt = Label(Cadre_4, text = 'Name of new file :', width=20)
file_name_2_value = StringVar() 
file_name_2 = Entry(Cadre_4, textvariable=file_name_2_value, width=20)

# bouton features
bouton_features=Button(Cadre_4, text="Get features ", command= lambda: ftrs.getFeaturesFromOrganismScript(file_name=file_name_2.get(), path=path_2.get(), organism_id=organism_id_3.get(), is_phage=is_phage_2_value.get()))

# display in the interface
Cadre_4.pack()
organism_id_3_txt.pack()
organism_id_3_txt.place(x=10, y=10)
organism_id_3.pack()
organism_id_3.place(x=150, y=10)
path_2_txt.pack()
path_2_txt.place(x=10, y=30)
path_2.pack()
path_2.place(x=150, y=30)
file_name_2_txt.pack()
file_name_2_txt.place(x=10, y=50)
file_name_2.pack()
file_name_2.place(x=150, y=50)
bouton_features.pack()
bouton_features.place(x=10, y=70)

# ================================================================
# Display Graph
# ================================================================

# Cadre 5
Cadre_5 = LabelFrame(fenetre, text="Display Graph", height=100, width=550, relief=GROOVE, labelanchor='n')

# liste graphiques
liste_graphique_value = StringVar(Cadre_5)
liste_graphique_value.set("Graphic choice")
liste_graphique = OptionMenu(Cadre_5, liste_graphique_value, "hierarchique", "heatmap", "clustermap")

# path
path_3_txt = Label(Cadre_5, text = 'File + Path :', width=20)
path_3_value = StringVar() 
path_3 = Entry(Cadre_5, textvariable=path_3_value, width=20)

# bouton graphic
bouton_graphic=Button(Cadre_5, text="Display graph", command= lambda: display_graphic.DisplayGraphicFromMatrix(path_3_value.get(), liste_graphique_value.get()))

# bouton open path 
bouton_open_path=Button(Cadre_5, text="Choose csv file", command= lambda: openpath(path_3_value))

# display in the interface
Cadre_5.pack() 
path_3_txt.pack()
path_3_txt.place(x=10, y=10)
path_3.pack()
path_3.place(x=150, y=10)
liste_graphique.pack()
liste_graphique.place(x=300, y=30)
bouton_open_path.pack()
bouton_open_path.place(x=150, y=30)
bouton_graphic.pack()
bouton_graphic.place(x=10,y=30)

fenetre.mainloop()
from tkinter import *
from tkinter.filedialog import askopenfilename

import network_chart_bacterium_script as ncbs
import network_chart_couples_script as nccs
import matrix_similarity_script as mss
import dendrogram

def openpath(e:Entry):
    path = askopenfilename(title="Ouvrir une image",filetypes=[('png files','.png'),('all files','.*')])
    e.set(path)

# fenetre principale
fenetre = Tk()
fenetre.geometry('500x500')
label = Label(fenetre, text="V0.1")

# liste graphiques
liste_graphique_value = StringVar(fenetre)
liste_graphique_value.set("Graphic choice")
liste_graphique = OptionMenu(fenetre, liste_graphique_value, "hierarchique", "heatmap", "clustermap")

# liste lysis_couple
liste_lysis_couple_value = StringVar(fenetre)
liste_lysis_couple_value.set("Lysis Type")
liste_lysis_couple = OptionMenu(fenetre, liste_lysis_couple_value, "CLEAR_LYSIS", "SEMI_CLEAR_LYSIS", "OPAQUE_LYSIS", "CLEAR_LYSIS_1E7PLUS", 
                                                                    "SEMI_CLEAR_LYSIS_1E7PLUS", "CLEAR_LYSIS_1E7MINUS", "SEMI_CLEAR_LYSIS_1E7MINUS",
                                                                    "ALL_CLEAR_LYSIS", "ALL_SEMI_CLEAR_LYSIS")

# organism_id
organism_id_txt = Label(fenetre, text = 'organism id :', width=20)
organism_id_value = IntVar() 
organism_id = Entry(fenetre, textvariable=organism_id_value, width=20)

# checkbutton
is_phage_value = IntVar()
is_phage = Checkbutton(fenetre, text="Phage ?", variable=is_phage_value)

# phage_1_id
phage_1_id_txt = Label(fenetre, text = 'phage 1 id :', width=20)
phage_1_id_value = IntVar() 
phage_1_id = Entry(fenetre, textvariable=phage_1_id_value, width=20)

# phage_2_id
phage_2_id_txt = Label(fenetre, text = 'phage 2 id :', width=20)
phage_2_id_value = IntVar() 
phage_2_id = Entry(fenetre, textvariable=phage_2_id_value, width=20)

# path
path_txt = Label(fenetre, text = 'path :', width=20)
path_value = StringVar() 
path = Entry(fenetre, textvariable=path_value, width=20)

# file_name
file_name_txt = Label(fenetre, text = 'file name :', width=20)
file_name_value = StringVar() 
file_name = Entry(fenetre, textvariable=file_name_value, width=20)

# bouton network_chart_script
bouton_network_chart_script=Button(fenetre, text="graphic : bacterium - phage", command= lambda: ncbs.networkChartOrganismScript(organism_id=organism_id.get(), is_phage=is_phage_value.get()))

# bouton network_chart_couple_script
bouton_network_chart_couple_script=Button(fenetre, text="graphic COUPLES: bacterium - phage", command= lambda: nccs.networkChartCouplesScript(parameter_type=liste_lysis_couple_value.get()))

# bouton matrice
bouton_matrix=Button(fenetre, text="create similitary matrix", command= lambda: mss.matrixSimilarityScript(file_name=file_name.get(), path=path.get(), organism_id=organism_id.get(), is_phage=is_phage_value.get()))

# bouton dendogram
bouton_dendogram=Button(fenetre, text="display graph", command= lambda: dendrogram.DisplayDendrogramFromMatrix(file_name.get(), liste_graphique_value.get()))

# bouton open path 
bouton_open_path=Button(fenetre, text="open path file", command= lambda: openpath(path_value))

phage_1_id_txt.pack()
phage_1_id_txt.place(x=10, y=25)
phage_1_id.pack()
phage_1_id.place(x=150, y=25)

phage_2_id_txt.pack()
phage_2_id_txt.place(x=10, y=50)
phage_2_id.pack()
phage_2_id.place(x=150, y=50)

organism_id_txt.pack()
organism_id_txt.place(x=10, y=75)
organism_id.pack()
organism_id.place(x=150, y=75)

path_txt.pack()
path_txt.place(x=10, y=100)
path.pack()
path.place(x=150, y=100)

file_name_txt.pack()
file_name_txt.place(x=10, y=125)
file_name.pack()
file_name.place(x=150, y=125)

bouton_network_chart_script.pack()
bouton_network_chart_script.place(x=10, y=150)
bouton_network_chart_couple_script.pack()
bouton_network_chart_couple_script.place(x=10, y=175)
bouton_matrix.pack()
bouton_matrix.place(x=10, y=200)
bouton_dendogram.pack()
bouton_dendogram.place(x=10,y=225)
bouton_open_path.pack()
bouton_open_path.place(x=320, y=120)
liste_graphique.pack()
liste_graphique.place(x=320, y=215)
liste_lysis_couple.pack()
liste_lysis_couple.place(x=320, y=170)
is_phage.pack()
is_phage.place(x=320, y=75)


fenetre.mainloop()
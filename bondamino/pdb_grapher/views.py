import os
import uuid
import Bio.PDB
import networkx as nx
from pyvis.network import Network
import numpy as np
import warnings
import random
from django.shortcuts import render, redirect, get_object_or_404
from django.urls import reverse
from django.contrib import messages
from django.conf import settings
from django.core.files.base import ContentFile
from .forms import PDBUploadForm, GraphOptionsForm
from .models import PDBFile, GeneratedGraph


AMINO_ACID_BONDS = { # Gekürzt für Lesbarkeit
    'ALA': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB')],
    'ARG': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'NE'), ('NE', 'CZ'), ('CZ', 'NH1'), ('CZ', 'NH2')],
    'ASN': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'ND2')],
    'ASP': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'OD1'), ('CG', 'OD2')],
    'CYS': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'SG')],
    'GLN': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'NE2')],
    'GLU': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'OE1'), ('CD', 'OE2')],
    'GLY': [('N', 'CA'), ('CA', 'C'), ('C', 'O')],
    'HIS': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'ND1'), ('ND1', 'CE1'), ('CE1', 'NE2'), ('NE2', 'CD2'), ('CD2', 'CG')], # Ring
    'ILE': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG1'), ('CG1', 'CD1'), ('CB', 'CG2')],
    'LEU': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CG', 'CD2')],
    'LYS': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'CE'), ('CE', 'NZ')],
    'MET': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'SD'), ('SD', 'CE')],
    'PHE': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CD1', 'CE1'), ('CE1', 'CZ'), ('CZ', 'CE2'), ('CE2', 'CD2'), ('CD2', 'CG')], # Ring
    'PRO': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD'), ('CD', 'N')], # Spezielle N-CD Bindung im Ring
    'SER': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'OG')],
    'THR': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'OG1'), ('CB', 'CG2')],
    'TRP': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CD1', 'NE1'), ('NE1', 'CE2'), ('CE2', 'CZ2'), ('CZ2', 'CH2'), ('CH2', 'CZ3'), ('CZ3', 'CE3'), ('CE3', 'CD2'), ('CD2', 'CG'), ('CE2', 'CD2')], # Indolring
    'TYR': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG'), ('CG', 'CD1'), ('CD1', 'CE1'), ('CE1', 'CZ'), ('CZ', 'OH'), ('CZ', 'CE2'), ('CE2', 'CD2'), ('CD2', 'CG')], # Ring
    'VAL': [('N', 'CA'), ('CA', 'C'), ('C', 'O'), ('CA', 'CB'), ('CB', 'CG1'), ('CB', 'CG2')],
}

_used_colors = set()

def get_distinct_color():
     global _used_colors

     _used_colors.clear()

     while True:
         color = "#{:06x}".format(random.randint(0, 0xFFFFFF))
         r, g, b = int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)
         brightness = (r * 299 + g * 587 + b * 114) / 1000
         if 50 < brightness < 200 and color not in _used_colors:
              _used_colors.add(color)
              return color
         # Fallback
         if len(_used_colors) > 500: return "#808080"


warnings.filterwarnings("ignore", category=Bio.PDB.PDBExceptions.PDBConstructionWarning)

def get_model_count(pdb_filepath):
    """Determines the number of models in a PDB file."""

    parser = Bio.PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("temp_structure", pdb_filepath)
        return len(structure)

    except Exception as e:
        print(f"Error while determining number of models in {pdb_filepath}: {e}")
        return 0

def create_protein_graph(pdb_filepath, model_id=0, label_type='atom_name', output_dir=None):
    """
    Reads PDB, creates a graph and saves the result as PyVis HTML in specified directory.
    Returns the path to the generated HTML file or None if an error occurs.
    """

    parser = Bio.PDB.PDBParser()
    try:
        structure = parser.get_structure("protein", pdb_filepath)
    except Exception as e:
        print(f"Error while parsing {pdb_filepath}: {e}")
        return None, f"Error while parsing PDB file: {e}"


    num_models = len(structure)

    if num_models == 0:
        return None, "Error: PDB file has no models."
    if model_id >= num_models:
        print(f"Warning: Model ID {model_id} is invalid for {num_models} models. Use model 0.")
        model_id = 0
    model = structure[model_id]

    # Initialize Graph
    G = nx.Graph()

    # generate unique file name for HTML output
    graph_filename = f"{uuid.uuid4().hex}.html"

    # check that directory exists
    os.makedirs(output_dir, exist_ok=True)
    output_html_path = os.path.join(output_dir, graph_filename)

    pv_network = Network(notebook=False, height='800px', width='100%', bgcolor='#ffffff', font_color='black', directed=False)

    # color assignments for nodes and edges
    residue_colors = {}
    unique_residues = set(residue.get_resname() for chain in model for residue in chain if Bio.PDB.is_aa(residue, standard=True))
    for res_name in unique_residues:
        residue_colors[res_name] = get_distinct_color()
    residue_colors['UNKNOWN'] = '#808080'

    # add nodes
    atom_dict = {}
    nodes_added = 0

    for chain in model:
        for residue in chain:
            # determine color by residue name
            res_name = residue.get_resname()
            if not Bio.PDB.is_aa(residue, standard=True) and res_name in AMINO_ACID_BONDS: standard_res_name = res_name
            elif Bio.PDB.is_aa(residue, standard=True): standard_res_name = res_name
            else: standard_res_name = "UNKNOWN"
            res_id = f"{chain.id}_{residue.get_id()[1]}"
            color = residue_colors.get(standard_res_name, residue_colors['UNKNOWN'])

            for atom in residue:
                # get atom details, node id, label and so on
                atom_name = atom.get_name().strip()
                element = atom.element.strip().upper() if atom.element else '?'
                coords = atom.get_coord()
                node_id = f"{res_id}_{atom_name}"
                node_label = atom_name if label_type == 'atom_name' else element
                title = (f"Atom: {atom_name}\n"
                         f"Element: {element}\n"
                         f"Residue: {res_name} {residue.get_id()[1]}\n"
                         f"Chain: {chain.id}\n"
                         f"Coords: ({coords[0]:.2f}, {coords[1]:.2f}, {coords[2]:.2f})")

                G.add_node(node_id,
                           label=node_label,
                           title=title,
                           color=color,
                           element=element,
                           atom_name=atom_name,
                           residue_name=res_name,
                           residue_id=res_id,
                           coords=coords)

                pv_network.add_node(node_id, label=node_label, title=title, color=color, shape='dot', size=10)
                atom_dict[node_id] = atom
                nodes_added += 1

    # add edges
    added_edges = set()
    edges_added = 0
    residue_list = list(model.get_residues())

    for i, residue in enumerate(residue_list):
        res_name = residue.get_resname()
        chain_id = residue.get_parent().id
        res_seq_id = residue.get_id()[1]
        current_res_id_str = f"{chain_id}_{res_seq_id}"

        # intra-residue bondings
        if res_name in AMINO_ACID_BONDS:
            bond_template = AMINO_ACID_BONDS[res_name]
            for atom1_name, atom2_name in bond_template:
                node1_id = f"{current_res_id_str}_{atom1_name}"
                node2_id = f"{current_res_id_str}_{atom2_name}"

                # check if both atoms exist in the model
                if node1_id in atom_dict and node2_id in atom_dict:
                    atom1 = atom_dict[node1_id]
                    atom2 = atom_dict[node2_id]

                    # avoid that the edge is not added twice
                    edge_tuple = tuple(sorted((node1_id, node2_id)))
                    if edge_tuple not in added_edges:
                        distance = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
                        edge_label = f"{distance:.2f} Å"

                        G.add_edge(node1_id, node2_id, length=distance, label=edge_label)
                        pv_network.add_edge(node1_id, node2_id, title=edge_label, value=distance,
                                             label=edge_label if len(G.edges) < 200 else None)
                        added_edges.add(edge_tuple)

        # peptide bonds (between C of this residue and N from the next residue)
        if Bio.PDB.is_aa(residue, standard=True) and i + 1 < len(residue_list):
            next_residue = residue_list[i+1]
            # check, if the next residue is in the same chain and that it is an amino acid
            if next_residue.get_parent().id == chain_id and Bio.PDB.is_aa(next_residue, standard=True):
                # check, if sequence numbers are Prüfen ob Sequenznummern follow one after the other
                if next_residue.get_id()[1] == res_seq_id + 1:
                    atom_C_name = 'C'
                    atom_N_name = 'N'
                    node_C_id = f"{current_res_id_str}_{atom_C_name}"
                    next_res_id_str = f"{chain_id}_{next_residue.get_id()[1]}"
                    node_N_id = f"{next_res_id_str}_{atom_N_name}"

                    if node_C_id in atom_dict and node_N_id in atom_dict:
                        atom_C = atom_dict[node_C_id]
                        atom_N = atom_dict[node_N_id]

                        edge_tuple = tuple(sorted((node_C_id, node_N_id)))
                        if edge_tuple not in added_edges:
                            distance = np.linalg.norm(atom_C.get_coord() - atom_N.get_coord())
                            if 1.2 < distance < 1.5:
                                edge_label = f"{distance:.2f} Å"
                                G.add_edge(node_C_id, node_N_id, length=distance, label=edge_label)
                                pv_network.add_edge(node_C_id, node_N_id, title=edge_label, value=distance,
                                                    # peptide bond color
                                                     label=edge_label if len(G.edges) < 200 else None, color = '#A0A0A0')
                                added_edges.add(edge_tuple)

        pass

    # adapt visualization and save it
    print(f"Graph for model {model_id} created: {nodes_added} nodes, {edges_added} edges")
    if nodes_added == 0:
        return None, f"No atoms in selected model ({model_id})."

    pv_network.set_options("""
        var options = {
          "nodes": {
            "font": {
              "size": 10
            }
          },
          "edges": {
            "color": {
              "inherit": true // Kantenfarbe von Knoten erben (kann überschrieben werden)
            },
            "smooth": { // Kann helfen bei überlappenden Kanten, kann aber langsam sein
               "enabled": true,
               "type": "dynamic", // "continuous" oder "dynamic"
               "roundness": 0.5
             },
             "font": {
                  "size": 8,
                  "align": "top" // Position des Kantenlabels
                }
          },
          "physics": { // Layout-Algorithmus
            "enabled": true, // Aktiviert Physik-Simulation für Layout
            "barnesHut": { // Ein gängiger Algorithmus
              "gravitationalConstant": -8000,
              "centralGravity": 0.3,
              "springLength": 100, // Beeinflusst nicht direkt die Kantenlänge basierend auf 'length'
              "springConstant": 0.04,
              "damping": 0.09,
              "avoidOverlap": 0.1 // Versucht Knotenüberlappung zu vermeiden
            },
            "solver": "barnesHut", // Alternativen: "forceAtlas2Based", "repulsion"
             "stabilization": { // Versucht das Layout zu stabilisieren
                "enabled": true,
                "iterations": 1000,
                "updateInterval": 50
              }
          },
          "interaction": {
            "tooltipDelay": 200, // Verzögerung für Tooltips (ms)
            "hideEdgesOnDrag": true, // Kanten beim Ziehen ausblenden (Performance)
            "navigationButtons": true, // Zoom/Navigationsbuttons
            "keyboard": true // Tastaturnavigation
          }
        }
        """)

    try:
        pv_network.save_graph(output_html_path)
        print(f"Interactive Graph temporarily saved as: {output_html_path}")
        return output_html_path, None
    except Exception as e:
        print(f"Error while saving HTML file {output_html_path}: {e}")
        # delete if there is and partially created file
        if os.path.exists(output_html_path):
             try: os.remove(output_html_path)
             except OSError: pass
        return None, f"Error while saving the graph HTML file: {e}"

def index_view(request):
    if request.method == 'POST':
        form = PDBUploadForm(request.POST, request.FILES)
        if form.is_valid():
            uploaded_file = request.FILES['pdb_file']
            # New PDBFile object gets created (not saved)
            pdb_instance = PDBFile(
                original_filename=uploaded_file.name,
                # id will be automatic created
            )
            # save file temporarily, to check model number
            pdb_instance.pdb_file.save(uploaded_file.name, uploaded_file, save=False) # save=False, da wir erst num_models brauchen

            try:
                num_models = get_model_count(pdb_instance.pdb_file.path)
                if num_models > 0:
                    pdb_instance.num_models = num_models
                    pdb_instance.save() # Jetzt PdbFile Objekt mit num_models speichern
                    messages.success(request, f'File "{pdb_instance.original_filename}" successfully uploaded.')
                    return redirect('pdb_grapher:select_options', pdb_file_id=pdb_instance.id)
                else:
                    messages.error(request, f'Could not find models in file "{uploaded_file.name}" or file is corrupted.')
                    # Important: Since pdb_instance has not yet been saved, we have to delete the temporarily uploaded file manually
                    if pdb_instance.pdb_file and os.path.exists(pdb_instance.pdb_file.path):
                       os.remove(pdb_instance.pdb_file.path)

            except Exception as e:
                 messages.error(request, f'Error when processing the PDB file: {e}')
                 if pdb_instance.pdb_file and os.path.exists(pdb_instance.pdb_file.path):
                    os.remove(pdb_instance.pdb_file.path)

        # If form invalid or error during parsing -> display form again
        # The errors from form.errors or messages are displayed in the template
    else: # GET Request
        form = PDBUploadForm()
    return render(request, 'pdb_grapher/index.html', {'form': form})


def select_options_view(request, pdb_file_id):
    """Shows the selection for model and label type."""
    pdb_instance = get_object_or_404(PDBFile, pk=pdb_file_id)

    if request.method == 'POST':
        # The form is sent to ‘process_pdb_view’
        # This view should not actually handle POST
        return redirect('pdb_grapher:select_options', pdb_file_id=pdb_instance.id)
    else:
        # GET Request
        # Create form with dynamic choices for model_id
        form = GraphOptionsForm(num_models=pdb_instance.num_models)

    return render(request, 'pdb_grapher/select_options.html', {
        'form': form,
        'pdb_instance': pdb_instance
    })


def process_pdb_view(request, pdb_file_id):
    """Processes the option selection and starts the graph generation."""
    pdb_instance = get_object_or_404(PDBFile, pk=pdb_file_id)

    if request.method == 'POST':
        # Form must be initialised again with num_models for validation
        form = GraphOptionsForm(request.POST, num_models=pdb_instance.num_models)
        if form.is_valid():
            model_id = int(form.cleaned_data['model_id'])
            label_type = form.cleaned_data['label_type']

            # Path to the temporary directory for the generated HTML
            temp_graph_dir = os.path.join(settings.MEDIA_ROOT, 'temp_graphs')

            # create Graph
            generated_html_path, error_message = create_protein_graph(
                pdb_filepath=pdb_instance.pdb_file.path,
                model_id=model_id,
                label_type=label_type,
                output_dir=temp_graph_dir
            )

            if generated_html_path:
                # Successful, create GeneratedGraph instance
                graph_instance = GeneratedGraph(
                    pdb_file=pdb_instance,
                    model_id_used=model_id,
                    label_type_used=label_type
                    # id wird automatisch generiert
                )
                # Open the temporarily generated file and save it in the FileField
                try:
                    with open(generated_html_path, 'rb') as f:
                        # Set file name for the FileField
                        graph_filename = os.path.basename(generated_html_path)
                        graph_instance.graph_html_file.save(graph_filename, ContentFile(f.read()), save=True)

                    messages.success(request, "Graph successfully generated.")
                    # Delete temporary file
                    os.remove(generated_html_path)
                    # Continue to the results page
                    return redirect('pdb_grapher:show_results', graph_id=graph_instance.id)

                except Exception as e:
                     messages.error(request, f"Error when saving the graph: {e}")
                     # Delete temporary file if necessary
                     if os.path.exists(generated_html_path): os.remove(generated_html_path)

            else:
                messages.error(request, f'Error during graph creation: {error_message}')

        # If form invalid or error -> Back to selection page
        # Form must be recreated for GET rendering
        return redirect('pdb_grapher:select_options', pdb_file_id=pdb_instance.id)

    else: # GET request not allowed for this view
        return redirect('pdb_grapher:select_options', pdb_file_id=pdb_instance.id)


def show_results_view(request, graph_id):
    """Shows the results page with the embedded graph."""
    graph_instance = get_object_or_404(GeneratedGraph, pk=graph_id)
    return render(request, 'pdb_grapher/results.html', {'graph': graph_instance})
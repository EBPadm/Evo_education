import streamlit as st
import random
import graphviz
import re

st.title("DNA塩基の進化を再現してみよう")

def edge_len_mod(id):
    name = st.session_state.seqnames[id]
    bodylist = st.session_state.graph.body
    pattern = re.compile(fr"-> {name} ")
    matched = [(i, value) for i, value in enumerate(bodylist) if pattern.search(value)]
    text = matched[0][1]
    lenmatch = re.search(r"len=(\d+\.?\d*)", text)
    number = float(lenmatch.group(1))+ 0.2
    new_text = re.sub(r"(len=)\d+\.?\d*", fr"\g<1>{number}", text)
    bodylist[matched[0][0]] = new_text
    return(bodylist)


def color_sequence(sequence, highlight_indices):
    colored_sequence = ""
    
    # Loop through the sequence and color specific letters
    for i, letter in enumerate(sequence):
        if i in highlight_indices:
            colored_sequence += f"<span style='color:red;font-weight:bold;'>{letter}</span>"
        else:
            colored_sequence += f"<span style='color:black;'>{letter}</span>"
    
    return colored_sequence

def single_mutation(id):
    seq = st.session_state.sequences[id]
    # Randomly choose a position to substitute
    random_position = random.randint(0, len(seq) - 1)
    st.session_state.nucleotides[id].append(random_position)
    # Get the current nucleotide at that index
    original_nucleotide = seq[random_position]
    # Create a list of possible nucleotides excluding the original one
    possible_nucleotides = ['A', 'T', 'G', 'C']
    possible_nucleotides.remove(original_nucleotide)
    # Randomly choose a new nucleotide from the remaining options
    new_nucleotide = random.choice(possible_nucleotides)
    # Substitute the letter at the chosen position
    st.session_state.sequences[id] = seq[:random_position] + new_nucleotide + seq[random_position+1:]
    edge_len_mod(id)

def all_mutation():
    for i, seq in enumerate(st.session_state.sequences):
        single_mutation(i)

@st.dialog('種数の上限(14)に達しました')
def error_dialog(): 
    if st.button('OK'):
        st.rerun()

def speciation(id):
    if(st.session_state.seqnum == 14):
        error_dialog()
    else:
        if(st.session_state.seqnum == 1):
            ancestor = "MRCA"
            offspring1 = "A"
            offspring2 = "B"
        else:
            last_id=len(st.session_state.seqnames)
            ancestor = st.session_state.seqnames[id]
            offspring1 = ancestor + alphabet_list[last_id*2-2]
            offspring2 = ancestor + alphabet_list[last_id*2-1]
        st.session_state.graph.edge(ancestor,offspring1,len="0")
        st.session_state.graph.edge(ancestor,offspring2,len="0")
        st.session_state.sequences.insert(id+1,st.session_state.sequences[id])
        last_id=len(st.session_state.seqnames)
        original_name = st.session_state.seqnames[id]
        st.session_state.seqnames[id] = original_name + alphabet_list[last_id*2-2]
        st.session_state.seqnames.insert(id+1,original_name + alphabet_list[last_id*2-1])
        st.session_state.nucleotides.insert(id+1,st.session_state.nucleotides[id].copy())
        st.session_state.seqnum += 1

def writing_fasta(output_space):
    outtext = ""
    for i in range(len(st.session_state.seqnames)):
        outtext += f">{st.session_state.seqnames[i]}\n"
        outtext += f"{st.session_state.sequences[i]}\n"
    output_space.text(outtext)

def game_process(seq_length):
    if "graph" not in st.session_state:
        st.session_state.graph = graphviz.Digraph(engine='neato')
    if "seqnum" not in st.session_state:
        st.session_state.seqnum = 1
    if "seqnames" not in st.session_state:
        st.session_state.seqnames = ['']
    if "sequences" not in st.session_state:
        st.session_state.sequences = [''.join(random.choices('ATGC', k=seq_length))]
    if "nucleotides" not in st.session_state:
        st.session_state.nucleotides = [[]]
    if "start" not in st.session_state:
        st.session_state.start = 0
    if st.session_state.start == 0:
        first_but = st.empty()
        first_seq = st.empty()
        first_seq.markdown(f"<p>{st.session_state.sequences[0]}</p>", unsafe_allow_html=True)
        if first_but.button("種分化!"):
            st.session_state.start = 1
            first_seq.empty()
            speciation(0)
            st.rerun()
    else:
        st.button("全部に突然変異",on_click=all_mutation)
        for index in range(st.session_state.seqnum):
            col1, col2, col3, _ = st.columns([2,2,2,4])
            col1.markdown(f"<p style='font-size: 20px;'>{st.session_state.seqnames[index]}</p>", unsafe_allow_html=True)
            col2.button("種分化"+str(index+1),on_click=speciation,args=[index])
            col3.button("突然変異"+str(index+1),on_click=single_mutation,args=[index])
            colored_seq = color_sequence(st.session_state.sequences[index], st.session_state.nucleotides[index])
            st.markdown(f"<p>{colored_seq}</p>", unsafe_allow_html=True)  # Display the sequence with color
        st.graphviz_chart(st.session_state.graph)
        fasta_but = st.empty()
        res_but = st.empty()
        output_space=st.empty()
        if fasta_but.button("FASTA出力"):
            writing_fasta(output_space)
        if res_but.button("リセット"):
            st.session_state.seqnum = 1
            st.session_state.seqnames = ['']
            st.session_state.sequences = [''.join(random.choices('ATGC', k=seq_length))]
            st.session_state.nucleotides = [[]]
            st.session_state.graph = graphviz.Digraph(engine='neato')

seq_length = 300
alphabet_list = [chr(i) for i in range(ord('A'), ord('Z') + 1)]
game_process(seq_length)

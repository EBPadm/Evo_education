import streamlit as st
import pandas as pd
import graphviz

def parse_fasta(fasta_string):
    seq_names = []
    sequences = []
    current_seq = []
    for line in fasta_string.strip().split("\n"):
        if line.startswith(">"):  # Header line
            if current_seq:  # Store the previous sequence
                sequences.append("".join(current_seq))
                current_seq = []
            seq_names.append(line[1:].strip())  # Store sequence name
        else:
            current_seq.append(line.strip())  # Collect sequence lines

    if current_seq:  # Store the last sequence
        sequences.append("".join(current_seq))
    return seq_names, sequences

# Function to calculate the number of differences between two sequences
def count_differences(seq1, seq2):
    # Ensure both sequences are the same length
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length")
    # Calculate the number of differences
    return sum(1 for a, b in zip(seq1, seq2) if a != b)

def qmatrix_generation(df):
    qmat = df
    nrow = len(df)
    for i in range(nrow):
        for j in range(nrow):
            if i == j:
                qmat.iat[i,j]="-"
            else:
                ilist = list(df.iloc[i,:i]) + list(df.iloc[i,i+1:])
                jlist = list(df.iloc[:j,j]) + list(df.iloc[j+1:,j])
                ilist = [int(s) for s in ilist]
                jlist = [int(s) for s in jlist]
                qmat.iat[i,j]=(nrow-2)*df.iat[i,j]-sum(ilist)-sum(jlist)
    return qmat

def s0calculation(df):
    nrow = len(df)
    ssum = 0
    for k in range(nrow-1):
        for l in range(k+1,nrow):
            ssum += df.iat[k,l]
    return ssum / (nrow - 1)

def scalculation(df,i,j):
    nrow = len(df)
    ilist = list(df.iloc[i,:i]) + list(df.iloc[i,i+1:])
    jlist = list(df.iloc[:j,j]) + list(df.iloc[j+1:,j])
    ilist = [int(s) for s in ilist]
    jlist = [int(s) for s in jlist]
    ssum = 0
    for k in range(nrow-1):
        for l in range(k+1,nrow):
            ssum += df.iat[k,l]
    nt = nrow - 2
    return (df.iat[i,j] - sum(ilist)/nt - sum(jlist)/nt + ssum*2/nt) / 2

def smatrix_generation(df):
    smat = df.copy()
    nrow = len(df)
    minis = scalculation(df,0,1)
    mini = 0
    minj = 1
    for i in range(nrow):
        for j in range(nrow):
            if i == j:
                smat.iat[i,j]="-"
            else:
                smat.iat[i,j]=scalculation(df,i,j)
                if smat.iat[i,j] < minis:
                    minis = smat.iat[i,j]
                    mini = i
                    minj = j
    return smat, mini, minj

def tree_generation(df,defrom,deto,delen,es):
    nrow = len(df)
    enum = len(defrom)
    graph = graphviz.Graph(engine='neato')
    ssum = 0
    for k in range(nrow-1):
        for l in range(k+1,nrow):
            ssum += df.iat[k,l]
    relen = ssum / nrow / (nrow - 1)
    for taxon in df.index:
        graph.edge("Node0",taxon,len=str(int(relen*es)))
    for k in range(enum):
        graph.edge(defrom[k],deto[k],len=str(int(delen[k]*es)))
    return graph

def df_update(df,mini,minj,num):
    nrow = len(df)
    newdf = df.copy()
    newdf = newdf.drop(index = df.index[mini])
    newdf = newdf.drop(index = df.index[minj])
    newdf = newdf.drop(columns=[df.columns[mini]])
    newdf = newdf.drop(columns=[df.columns[minj]])
    newname = "Node" + str(int(num))
    newnode_d = []
    for k in range(nrow):
        if k != mini and k!= minj:
            newnode_d.append((df.iat[mini,k]+df.iat[minj,k]-df.iat[mini,minj])/2)
    newdf[newname]=newnode_d
    newnode_d.append("-")
    newdf.loc[newname]=newnode_d
    return newdf

def elencalculation(df,i,j):
    nrow = len(df)
    ilist = list(df.iloc[i,:i]) + list(df.iloc[i,i+1:])
    jlist = list(df.iloc[:j,j]) + list(df.iloc[j+1:,j])
    ilist = [int(s) for s in ilist]
    jlist = [int(s) for s in jlist]
    eleni = (df.iat[i,j] + (sum(ilist) - sum(jlist)) / (nrow - 2))/2
    elenj = df.iat[i,j] - eleni
    return eleni, elenj

escale = 0.5
st.title("NJ法で系統樹を作成する")
if "analysis" not in st.session_state:
    st.session_state.analysis = 0
if "process" not in st.session_state:
    st.session_state.process = 0
st.session_state.fasta_string = st.text_area("DNA配列の入力(fasta形式)",
""">Exs1
ATGCGGCTCAATTG
>Exs2
ATGCGACGTACTTG
>Exs3
ACGCGACGCACTAG
>Exs4
ATGCCGCGGAATGG
>Exs5
ATGCCGCGGTATGG"""
)
if st.button("解析") or st.session_state.analysis == 1:
    st.session_state.analysis = 1
    seqnames, sequences = parse_fasta(st.session_state.fasta_string)
    if len(sequences)>1:
        # Prepare a list to store the differences
        data = []
        # Compare each pair of sequences
        for i in range(len(sequences)):
            data.append([])
            for j in range(len(sequences)):
                if i==j:
                    data[i].append("-")
                else:
                    diff = count_differences(sequences[i], sequences[j])
                    data[i].append(diff)
        if "df" not in st.session_state:
            st.session_state.df = pd.DataFrame(data, columns=seqnames)
            st.session_state.df.index = seqnames
        if "predf" not in st.session_state:
            st.session_state.predf = st.session_state.df.copy()
        next_button = st.empty()
        st.write("距離行列")
        st.dataframe(st.session_state.predf)
        if "det_efrom" not in st.session_state:
            st.session_state.det_efrom = []
        if "det_eto" not in st.session_state:
            st.session_state.det_eto = []
        if "det_elen" not in st.session_state:
            st.session_state.det_elen = []
        if "tree" not in st.session_state:
            st.session_state.tree = tree_generation(st.session_state.df,st.session_state.det_efrom,st.session_state.det_eto,st.session_state.det_elen,escale)
        if "svalue" not in st.session_state:
            st.session_state.svalue = s0calculation(st.session_state.df)
        escale = st.slider("枝の長さスケール", 0.05, 1.0, value=0.5)
        if st.button("Refresh"):
            st.session_state.tree = tree_generation(st.session_state.df,st.session_state.det_efrom,st.session_state.det_eto,st.session_state.det_elen,escale)
            st.rerun()
        st.graphviz_chart(st.session_state.tree)
        st.write(f"全ての枝の長さ, S = {st.session_state.svalue}")
        smatrix, mini, minj = smatrix_generation(st.session_state.predf)
        st.write("S-matrix")
        st.dataframe(smatrix)
        if len(st.session_state.predf) == 3:
            next_button.empty()
        else:
            if next_button.button("次へ"):
                st.session_state.process += 1
                if st.session_state.process % 2 == 1:
                    st.session_state.svalue = smatrix.iat[mini, minj]
                    node_num = (st.session_state.process - 1) / 2 + 1
                    st.session_state.det_efrom.append(f"Node{int(node_num)}")
                    st.session_state.det_eto.append(st.session_state.df.index[mini])
                    eleni, elenj = elencalculation(st.session_state.df,mini,minj)
                    st.session_state.det_elen.append(eleni)
                    st.session_state.det_efrom.append(f"Node{int(node_num)}")
                    st.session_state.det_eto.append(st.session_state.df.index[minj])
                    st.session_state.det_elen.append(elenj)
                    st.session_state.df = df_update(st.session_state.df,mini,minj,node_num)
                    st.session_state.tree = tree_generation(st.session_state.df,st.session_state.det_efrom,st.session_state.det_eto,st.session_state.det_elen,escale)
                else:
                    st.session_state.predf = st.session_state.df.copy()

            

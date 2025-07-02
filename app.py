import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from io import BytesIO
from Bio import SeqIO
import json
import plotly.express as px
import plotly.graph_objects as go
import seaborn as sns
import numpy as np
from matplotlib.patches import Rectangle
import math
import zipfile
import random
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import io
import openpyxl

# --- Add simple step-by-step instructions at the top of the app ---
def show_simple_steps():
    st.markdown("""
    ## How to Use This App (Step-by-Step)
    1. **Choose how to input your sequences** (CSV/Excel, FASTA, GenBank, or Manual Entry).
    2. **Upload your file or enter your sequences** in the main panel.
    3. **Click 'Calculate GC Content'** to analyze your sequences.
    4. **(Optional) Use the Sequence Fragmentation Tool** in the sidebar to split your sequences into smaller fragments for further analysis.
    5. **Download your results** in your preferred format (Excel, CSV, FASTA, GenBank, etc.).
    
    > **Tip:** You must upload or enter sequences before you can use the fragmentation tool!
    """)

# --- Sequence Sanitization and Validation ---
def sanitize_sequence(seq):
    if not isinstance(seq, str):
        seq = str(seq)
    return ''.join(filter(lambda x: x.upper() in ['A', 'T', 'G', 'C'], seq.upper()))

def is_valid_sequence(seq):
    if not isinstance(seq, str):
        seq = str(seq)
    return all(base in ['A', 'T', 'G', 'C'] for base in seq.upper())

# --- Sequence Fragmentation Tool ---
def fragment_sequence(name, seq, fragment_length):
    if not isinstance(seq, str):
        seq = str(seq)
    seq = sanitize_sequence(seq)
    if len(seq) == 0:
        return []
    fragments = []
    for i in range(0, len(seq), fragment_length):
        fragment = seq[i:i + fragment_length]
        if len(fragment) >= fragment_length:
            fragment_name = f"{name}_fragment_{i//fragment_length + 1}"
            fragments.append((fragment_name, fragment))
    return fragments

def process_sequences_for_fragmentation(sequences, fragment_length):
    all_fragments = []
    for name, seq in sequences:
        fragments = fragment_sequence(name, seq, fragment_length)
        all_fragments.extend(fragments)
    return all_fragments

# --- Nucleotide Analysis ---
def analyze_sequence(name, seq):
    if not isinstance(seq, str):
        seq = str(seq)
    seq = sanitize_sequence(seq)
    length = len(seq)
    if length == 0:
        return {"Gene Name": name, "Fragment Length": "Invalid sequence or empty after sanitization"}
    a_count = seq.count('A')
    t_count = seq.count('T')
    g_count = seq.count('G')
    c_count = seq.count('C')
    gc_count = g_count + c_count
    at_count = a_count + t_count
    return {
        "Gene Name": name,
        "Sequence": seq,
        "Length": length,
        "A Count": a_count,
        "T Count": t_count,
        "G Count": g_count,
        "C Count": c_count,
        "A %": round((a_count / length) * 100, 2),
        "T %": round((t_count / length) * 100, 2),
        "G %": round((g_count / length) * 100, 2),
        "C %": round((c_count / length) * 100, 2),
        "GC %": round((gc_count / length) * 100, 2),
        "AT %": round((at_count / length) * 100, 2),
    }

def process_fasta(file, max_sequences=1000):
    sequences = []
    try:
        file.seek(0)
        file_content = file.read()
        encodings = ['utf-8', 'latin-1', 'cp1252']
        text = None
        for encoding in encodings:
            try:
                text = file_content.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        if text is None:
            st.error("Could not decode the FASTA file. Please ensure it's a valid text file.")
            return []
        lines = text.split('\n')
        current_name = None
        current_sequence = []
        count = 0
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name and current_sequence:
                    full_sequence = ''.join(current_sequence)
                    if len(full_sequence) > 0:
                        sequences.append((current_name, full_sequence))
                        count += 1
                        if count >= max_sequences:
                            break
                current_name = line[1:].strip()
                current_sequence = []
            else:
                if current_name:
                    current_sequence.append(line)
        if current_name and current_sequence and count < max_sequences:
            full_sequence = ''.join(current_sequence)
            if len(full_sequence) > 0:
                sequences.append((current_name, full_sequence))
                count += 1
        if count == 0:
            st.warning("No sequences found in the uploaded FASTA file. Please check the file format.")
        elif count >= max_sequences:
            st.info(f"Processed the first {max_sequences} sequences from the FASTA file.")
        else:
            st.success(f"Successfully processed {count} sequences from the FASTA file.")
    except Exception as e:
        st.error(f"Error parsing FASTA file: {str(e)}")
        st.info("Please ensure your FASTA file has the correct format with '>' headers followed by sequence data.")
        return []
    return sequences

def process_genbank(file, max_sequences=1000):
    sequences = []
    try:
        file.seek(0)
        file_content = file.read()
        encodings = ['utf-8', 'latin-1', 'cp1252']
        text = None
        for encoding in encodings:
            try:
                text = file_content.decode(encoding)
                break
            except UnicodeDecodeError:
                continue
        if text is None:
            st.error("Could not decode the GenBank file. Please ensure it's a valid text file.")
            return []
        text_io = io.StringIO(text)
        records = SeqIO.parse(text_io, "genbank")
        count = 0
        for record in records:
            if count >= max_sequences:
                break
            seq_str = str(record.seq)
            if len(seq_str) > 0:
                sequences.append((str(record.id), seq_str))
                count += 1
        if count == 0:
            st.warning("No sequences found in the uploaded GenBank file. Please check the file format.")
        elif count >= max_sequences:
            st.info(f"Processed the first {max_sequences} sequences from the GenBank file.")
        else:
            st.success(f"Successfully processed {count} sequences from the GenBank file.")
    except Exception as e:
        st.error(f"Error parsing GenBank file: {str(e)}")
        st.info("Please ensure your GenBank file has the correct format with LOCUS, DEFINITION, and ORIGIN sections.")
        return []
    return sequences

# ... (rest of your app logic, including export functions, main(), etc.) ...
# Please let me know if you want the full code pasted here, or if you want to restore from a backup file.

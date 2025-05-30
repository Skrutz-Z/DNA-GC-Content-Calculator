import streamlit as st
import pandas as pd
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

# --- Instructions ---
def show_instructions():
    st.markdown("""
    ### 📋 Instructions
    
    1. **Choose Input Method**:
       - **Upload CSV/Excel**: Upload a CSV or Excel file with 'Gene Name' and 'Sequence' columns
       - **Upload FASTA**: Upload a FASTA file containing DNA sequences
       - **Manual Entry**: Enter sequences manually with gene names
    
    2. **Input Format**:
       - Only A, T, G, C nucleotides are accepted (case insensitive)
       - Invalid characters will be automatically removed
       - Each sequence must have a unique gene name
    
    3. **Analysis**:
       - Click "Calculate GC Content" to process sequences
       - View results in the interactive table
       - Explore visualizations of GC content and nucleotide composition
    
    4. **Export Options**:
       - Export results in multiple formats (Excel, CSV, JSON)
       - Use "Export All Results" for a complete dataset export
       - Customize output filename before downloading
    """)

# --- Sequence Sanitization and Validation ---
def sanitize_sequence(seq):
    # Convert to string if not already
    if not isinstance(seq, str):
        seq = str(seq)
    return ''.join(filter(lambda x: x.upper() in ['A', 'T', 'G', 'C'], seq.upper()))

def is_valid_sequence(seq):
    # Convert to string if not already
    if not isinstance(seq, str):
        seq = str(seq)
    return all(base in ['A', 'T', 'G', 'C'] for base in seq.upper())

# --- Nucleotide Analysis ---
def analyze_sequence(name, seq):
    # Convert sequence to string if it's not already
    if not isinstance(seq, str):
        seq = str(seq)
    
    seq = sanitize_sequence(seq)
    length = len(seq)
    if length == 0:
        return {"Gene Name": name, "Error": "Invalid sequence or empty after sanitization"}

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

def process_fasta(file):
    sequences = []
    records = SeqIO.parse(file, "fasta")
    for record in records:
        sequences.append((record.id, str(record.seq)))
    return sequences

def display_visuals(df):
    st.subheader("📊 GC Content Distribution")
    fig, ax = plt.subplots()
    df.plot.bar(x='Gene Name', y='GC %', ax=ax, color='purple', legend=False)
    plt.ylabel('GC %')
    plt.xticks(rotation=45, ha='right')
    st.pyplot(fig)

    st.subheader("🧬 Nucleotide Composition Per Sequence")
    for _, row in df.iterrows():
        fig, ax = plt.subplots()
        ax.pie(
            [row['A %'], row['T %'], row['G %'], row['C %']],
            labels=['A %', 'T %', 'G %', 'C %'],
            autopct='%1.1f%%',
            startangle=90
        )
        ax.set_title(f"{row['Gene Name']} - Base % Composition")
        st.pyplot(fig)

def export_data(df, format_type):
    if format_type == "Excel":
        towrite = BytesIO()
        df.to_excel(towrite, index=False, engine='openpyxl')
        towrite.seek(0)
        return towrite, "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx"
    elif format_type == "CSV":
        towrite = BytesIO()
        df.to_csv(towrite, index=False)
        towrite.seek(0)
        return towrite, "text/csv", "csv"
    else:  # JSON
        towrite = BytesIO()
        json_str = df.to_json(orient='records', indent=2)
        towrite.write(json_str.encode())
        towrite.seek(0)
        return towrite, "application/json", "json"

def create_gc_heatmap(df):
    st.subheader("🌡️ GC Content Heatmap")
    
    # Create a matrix of GC content for each position
    sequences = df['Sequence'].tolist()
    max_len = max(len(seq) for seq in sequences)
    
    # Initialize matrix
    gc_matrix = np.zeros((len(sequences), max_len))
    
    # Fill matrix with GC content for each position
    for i, seq in enumerate(sequences):
        for j in range(len(seq)):
            if j < len(seq):
                window = seq[max(0, j-10):min(len(seq), j+11)]
                gc_count = window.count('G') + window.count('C')
                gc_matrix[i, j] = (gc_count / len(window)) * 100
    
    # Create heatmap using plotly
    fig = go.Figure(data=go.Heatmap(
        z=gc_matrix,
        x=list(range(max_len)),
        y=df['Gene Name'].tolist(),
        colorscale='Viridis',
        colorbar=dict(title='GC %')
    ))
    
    fig.update_layout(
        title='GC Content Distribution Across Sequences',
        xaxis_title='Position',
        yaxis_title='Gene Name',
        height=400 + (len(sequences) * 20)  # Adjust height based on number of sequences
    )
    
    st.plotly_chart(fig, use_container_width=True)

def calculate_information_content(freq):
    """Calculate information content in bits."""
    if freq == 0:
        return 0
    return freq * math.log2(freq * 4)  # 4 for number of nucleotides

def create_sequence_logo(sequences):
    st.subheader("🎨 Sequence Logo")
    
    # Calculate position frequency matrix
    max_len = max(len(seq) for seq in sequences)
    pfm = np.zeros((4, max_len))  # 4 nucleotides
    
    for seq in sequences:
        for i, base in enumerate(seq):
            if i < max_len:
                if base == 'A':
                    pfm[0, i] += 1
                elif base == 'T':
                    pfm[1, i] += 1
                elif base == 'G':
                    pfm[2, i] += 1
                elif base == 'C':
                    pfm[3, i] += 1
    
    # Normalize
    pfm = pfm / len(sequences)
    
    # Calculate information content
    ic = np.zeros(max_len)
    for i in range(max_len):
        ic[i] = sum(calculate_information_content(freq) for freq in pfm[:, i])
    
    # Create the plot
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Colors for nucleotides
    colors = {'A': '#2ecc71', 'T': '#e74c3c', 'G': '#f1c40f', 'C': '#3498db'}
    bases = ['A', 'T', 'G', 'C']
    
    # Plot each position
    for i in range(max_len):
        # Sort frequencies for this position
        freqs = pfm[:, i]
        sorted_indices = np.argsort(freqs)
        
        # Plot each base
        y_bottom = 0
        for idx in sorted_indices:
            if freqs[idx] > 0:
                height = freqs[idx] * ic[i]
                rect = Rectangle((i, y_bottom), 1, height,
                               facecolor=colors[bases[idx]],
                               edgecolor='black',
                               linewidth=0.5)
                ax.add_patch(rect)
                y_bottom += height
    
    # Customize the plot
    ax.set_xlim(0, max_len)
    ax.set_ylim(0, max(ic) * 1.1)
    ax.set_xlabel('Position')
    ax.set_ylabel('Bits')
    ax.set_title('Sequence Logo')
    
    # Add legend
    legend_elements = [Rectangle((0, 0), 1, 1, facecolor=color, edgecolor='black')
                      for color in colors.values()]
    ax.legend(legend_elements, bases, loc='upper right')
    
    # Remove spines
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    
    # Show the plot
    st.pyplot(fig)
    plt.close(fig)

def create_interactive_plots(df):
    st.subheader("📊 Interactive Plots")
    
    # GC Content Distribution
    fig_gc = px.box(df, y='GC %', title='GC Content Distribution',
                    hover_data=['Gene Name', 'Length', 'GC %'])
    fig_gc.update_traces(marker_color='purple')
    st.plotly_chart(fig_gc, use_container_width=True)
    
    # Nucleotide Composition
    fig_comp = px.bar(df, 
                     x='Gene Name',
                     y=['A %', 'T %', 'G %', 'C %'],
                     title='Nucleotide Composition by Gene',
                     barmode='group',
                     hover_data=['Length', 'GC %'])
    st.plotly_chart(fig_comp, use_container_width=True)
    
    # Length vs GC Content
    fig_scatter = px.scatter(df,
                           x='Length',
                           y='GC %',
                           color='GC %',
                           hover_data=['Gene Name', 'Length', 'GC %'],
                           title='Sequence Length vs GC Content')
    st.plotly_chart(fig_scatter, use_container_width=True)

def export_all_to_zip(df, output_filename):
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, "w") as zip_file:
        # Excel
        excel_buffer = BytesIO()
        df.to_excel(excel_buffer, index=False, engine='openpyxl')
        excel_buffer.seek(0)
        zip_file.writestr(f"{output_filename}.xlsx", excel_buffer.read())
        # CSV
        csv_buffer = BytesIO()
        df.to_csv(csv_buffer, index=False)
        csv_buffer.seek(0)
        zip_file.writestr(f"{output_filename}.csv", csv_buffer.read())
        # JSON
        json_str = df.to_json(orient='records', indent=2)
        zip_file.writestr(f"{output_filename}.json", json_str)
    zip_buffer.seek(0)
    return zip_buffer

def main():
    st.set_page_config(page_title="GC Content Calculator", page_icon=None, layout="wide")
    st.markdown("""
        <h1 style='text-align: center; color: #6c63ff;'>GC Content Calculator</h1>
    """, unsafe_allow_html=True)

    # Sidebar with logo, navigation, and about info
    st.sidebar.image("https://cdn-icons-png.flaticon.com/512/616/616494.png", width=100)
    st.sidebar.title("Navigation")
    with st.sidebar.expander("About this app"):
        st.write("""
        This tool calculates GC content and provides interactive visualizations for DNA sequences. 
        Upload your data or enter sequences manually to explore nucleotide composition, GC content, and more. 
        
        **Made by Shubh Rakesh Nahar, Troy University.**
        """)
    # Fun facts about DNA/genes/sequences
    facts = [
        "The human genome contains about 3 billion base pairs.",
        "GC content can affect the stability of DNA.",
        "Some bacteria have extremely high or low GC content.",
        "DNA was first isolated by Friedrich Miescher in 1869.",
        "GC-rich regions are often found near gene promoters.",
        "Genes are segments of DNA that code for proteins.",
        "The longest human gene is over 2.4 million base pairs long!",
        "Mitochondrial DNA is inherited only from your mother.",
        "Some viruses use RNA instead of DNA as their genetic material.",
        "The fruit fly has about 15,000 genes, while humans have about 20,000-25,000.",
        "DNA stands for Deoxyribonucleic Acid.",
        "The double helix structure of DNA was discovered in 1953.",
        "Some plants have much more DNA than humans!",
        "The GC content of a genome can be used to identify species.",
        "DNA can be extracted from almost any living thing, even ancient fossils!"
    ]
    if 'fun_fact_idx' not in st.session_state:
        st.session_state['fun_fact_idx'] = random.randint(0, len(facts)-1)
    if st.sidebar.button("Show another fun fact"):
        prev_idx = st.session_state['fun_fact_idx']
        new_idx = prev_idx
        while new_idx == prev_idx:
            new_idx = random.randint(0, len(facts)-1)
        st.session_state['fun_fact_idx'] = new_idx
        st.rerun()
    st.sidebar.success(f"Fun Fact: {facts[st.session_state['fun_fact_idx']]}")

    st.markdown("---")
    st.markdown("#### Upload your data or enter sequences manually below.")

    input_method = st.radio("Choose input method", ["Upload CSV/Excel", "Upload FASTA", "Manual Entry"])
    sequences = []

    if input_method == "Upload CSV/Excel":
        uploaded_file = st.file_uploader("Upload a CSV or Excel file with 'Gene Name' and 'Sequence' columns", 
                                       type=["csv", "xlsx", "xls"])
        if uploaded_file:
            try:
                if uploaded_file.name.endswith((".xlsx", ".xls")):
                    df = pd.read_excel(uploaded_file)
                else:
                    df = pd.read_csv(uploaded_file)
                if "Gene Name" in df.columns and "Sequence" in df.columns:
                    df['Sequence'] = df['Sequence'].fillna('').astype(str)
                    df['Gene Name'] = df['Gene Name'].fillna('').astype(str)
                    sequences = list(zip(df["Gene Name"], df["Sequence"]))
                else:
                    st.toast("File must contain 'Gene Name' and 'Sequence' columns.", icon=None, duration=2)
            except Exception as e:
                st.toast(f"Error reading file: {str(e)}", icon=None, duration=2)

    elif input_method == "Upload FASTA":
        fasta_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])
        if fasta_file:
            sequences = process_fasta(fasta_file)

    elif input_method == "Manual Entry":
        num = st.number_input("How many sequences would you like to enter?", min_value=1, max_value=50, value=1)
        manual_entries = []
        for i in range(num):
            st.markdown(f"**Sequence {i+1}**")
            name = st.text_input(f"Gene Name {i+1}", key=f"name_{i}")
            seq = st.text_area(f"Sequence {i+1}", key=f"seq_{i}")
            if name and seq:
                manual_entries.append((name, seq))
        sequences.extend(manual_entries)

    if sequences:
        if st.button("Calculate GC Content"):
            results = [analyze_sequence(name, seq) for name, seq in sequences]
            result_df = pd.DataFrame(results)

            if "Error" in result_df.columns:
                st.toast("Some sequences were invalid and skipped.", icon=None, duration=2)
                result_df = result_df.dropna(subset=["Length"])

            st.toast("Analysis complete!", icon=None, duration=2)
            st.markdown("---")
            st.markdown("#### Results & Visualizations")
            with st.expander("Show Data Table", expanded=True):
                st.dataframe(result_df, use_container_width=True)
            with st.expander("Show GC Content Heatmap"):
                create_gc_heatmap(result_df)
            with st.expander("Show Sequence Logo"):
                create_sequence_logo([seq for _, seq in sequences])
            with st.expander("Show Interactive Plots"):
                create_interactive_plots(result_df)
            st.markdown("---")
            st.subheader("Export Results")
            output_filename = st.text_input("Enter output file name (without extension):", "gc_output_v4")
            if st.button("Export All Results"):
                zip_buffer = export_all_to_zip(result_df, output_filename)
                st.toast("All results exported as ZIP!", icon=None, duration=2)
                st.download_button(
                    label="Download All Results (ZIP)",
                    data=zip_buffer,
                    file_name=f"{output_filename}_all_results.zip",
                    mime="application/zip"
                )
            excel_buffer = BytesIO()
            result_df.to_excel(excel_buffer, index=False, engine='openpyxl')
            excel_buffer.seek(0)
            st.download_button(
                label="Download Excel Only",
                data=excel_buffer,
                file_name=f"{output_filename}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
    st.markdown("---")
    st.markdown("<div style='text-align: center; color: #888;'>Made by Shubh Rakesh Nahar, Troy University </div>", unsafe_allow_html=True)

if __name__ == "__main__":
    main()

"""
=============================================================================
Bioinformatics II - Lab 1: Retrieval and Basic Analysis of Protein Sequence
Theme: Dark Biotech Interface (Glassmorphism & High-Contrast Neon Cyan)
=============================================================================
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Set Streamlit Page Configuration
st.set_page_config(
    page_title="Lab 1 - Protein Analysis",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ---------------------------------------------------------------------------
# CUSTOM CSS: DARK BIOTECH THEME
# ---------------------------------------------------------------------------
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500;700&display=swap');

    /* 1. Hide default Streamlit header and footer */
    #MainMenu {visibility: hidden;}
    header {visibility: hidden;}
    footer {visibility: hidden;}
    .stDeployButton {display:none;}

    /* 2. Global Dark Biotech Theme */
    .stApp {
        background: radial-gradient(circle at 15% 20%, #0d1527 0%, #070a13 55%, #04060a 100%) !important;
        color: #F8FAFC !important;
        font-family: 'Inter', sans-serif;
    }

    /* 3. Sidebar: Darker background with subtle neon cyan glowing border */
    section[data-testid="stSidebar"] {
        background-color: #050811 !important;
        border-right: 1px solid rgba(0, 255, 204, 0.28) !important;
        box-shadow: 4px 0 25px rgba(0, 255, 204, 0.08);
    }
    section[data-testid="stSidebar"] .stRadio > label {
        color: #00FFCC !important;
        font-weight: 700;
        font-size: 1rem;
        letter-spacing: 0.5px;
    }
    section[data-testid="stSidebar"] div[role="radiogroup"] label p {
        color: #FFFFFF !important;
        font-weight: 600;
    }

    /* 4. Glassmorphism Card Containers (Used for record info) */
    .glass-card {
        background: rgba(13, 22, 38, 0.75);
        backdrop-filter: blur(14px);
        -webkit-backdrop-filter: blur(14px);
        border: 1px solid rgba(0, 255, 204, 0.28);
        border-radius: 12px;
        padding: 22px;
        box-shadow: 0 8px 32px 0 rgba(0, 0, 0, 0.5), inset 0 0 12px rgba(0, 255, 204, 0.04);
        margin-bottom: 22px;
    }

    /* 5. Typography */
    .biotech-title {
        font-size: 2.2rem;
        font-weight: 700;
        color: #FFFFFF !important;
        letter-spacing: -0.5px;
        margin-bottom: 4px;
        text-shadow: 0 0 16px rgba(0, 255, 204, 0.35);
    }
    .biotech-subtitle {
        font-size: 1rem;
        color: #E2E8F0 !important;
        margin-bottom: 24px;
        letter-spacing: 0.2px;
    }
    .section-heading {
        color: #00FFCC !important;
        font-size: 1.25rem;
        font-weight: 700;
        letter-spacing: 0.5px;
        margin-top: 14px;
        margin-bottom: 14px;
        text-shadow: 0 0 10px rgba(0, 255, 204, 0.3);
    }
    .field-label {
        color: #E2E8F0 !important;
        font-size: 0.85rem;
        font-weight: 700;
        text-transform: uppercase;
        letter-spacing: 0.8px;
        display: block;
        margin-bottom: 4px;
    }

    /* 6. "Enter Uniprot ID" Font Colour: PURE WHITE */
    div[data-testid="stTextInput"] label p,
    div[data-testid="stTextInput"] label {
        color: #FFFFFF !important;
        font-size: 1.05rem !important;
        font-weight: 700 !important;
        letter-spacing: 0.5px !important;
    }

    /* 7. Text Input Box Styling */
    div[data-baseweb="input"] {
        background-color: rgba(9, 14, 26, 0.9) !important;
        border: 1px solid rgba(0, 255, 204, 0.35) !important;
        border-radius: 6px !important;
    }
    div[data-baseweb="input"]:focus-within {
        border-color: #00FFCC !important;
        box-shadow: 0 0 14px rgba(0, 255, 204, 0.45) !important;
    }
    div[data-baseweb="input"] input {
        color: #FFFFFF !important;
        font-size: 1.05rem !important;
        font-family: 'JetBrains Mono', monospace !important;
        font-weight: 600 !important;
    }

    /* 8. Neon Primary Button */
    .stButton>button {
        background: linear-gradient(135deg, rgba(0, 255, 204, 0.22), rgba(0, 255, 204, 0.08)) !important;
        color: #00FFCC !important;
        border: 1px solid #00FFCC !important;
        border-radius: 6px !important;
        font-weight: 700 !important;
        font-size: 0.95rem !important;
        letter-spacing: 0.5px !important;
        box-shadow: 0 0 16px rgba(0, 255, 204, 0.22) !important;
        transition: all 0.25s ease-in-out !important;
        margin-top: 4px;
    }
    .stButton>button:hover {
        background: #00FFCC !important;
        color: #050811 !important;
        box-shadow: 0 0 28px rgba(0, 255, 204, 0.65) !important;
        border-color: #00FFCC !important;
        transform: translateY(-1px);
    }

    /* 9. High-Contrast Prominent Stat Cards */
    .stat-card {
        background: linear-gradient(135deg, rgba(16, 28, 48, 0.85), rgba(9, 16, 30, 0.9));
        border: 1px solid rgba(0, 255, 204, 0.38);
        border-radius: 12px;
        padding: 18px 16px;
        text-align: center;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.4), inset 0 0 12px rgba(0, 255, 204, 0.05);
        transition: transform 0.2s ease, border-color 0.2s ease, box-shadow 0.2s ease;
    }
    .stat-card:hover {
        border-color: #00FFCC;
        transform: translateY(-2px);
        box-shadow: 0 0 24px rgba(0, 255, 204, 0.25);
    }
    .stat-label {
        color: #FFFFFF !important;
        font-size: 0.88rem !important;
        font-weight: 700 !important;
        text-transform: uppercase !important;
        letter-spacing: 1px !important;
        margin-bottom: 8px !important;
    }
    .stat-value {
        color: #00FFCC !important;
        font-family: 'JetBrains Mono', monospace !important;
        font-size: 2.1rem !important;
        font-weight: 700 !important;
        text-shadow: 0 0 16px rgba(0, 255, 204, 0.45) !important;
        line-height: 1.2 !important;
    }
    .stat-sub {
        color: #CBD5E1 !important;
        font-size: 0.82rem !important;
        margin-top: 6px !important;
        font-weight: 500 !important;
    }

    /* 10. FASTA Code Box Line-Wrapping (Full View Without Horizontal Scroll) */
    div[data-testid="stExpander"] pre,
    div[data-testid="stExpander"] code,
    pre, code {
        font-family: 'JetBrains Mono', 'Fira Code', 'Courier New', monospace !important;
        background-color: #070B14 !important;
        color: #38BDF8 !important;
        border: 1px solid rgba(56, 189, 248, 0.32) !important;
        border-radius: 6px !important;
        white-space: pre-wrap !important;
        word-break: break-all !important;
        overflow-x: hidden !important;
        line-height: 1.65 !important;
        font-size: 0.95rem !important;
    }

    /* Dropdown Expander Styling */
    details {
        background: rgba(13, 22, 38, 0.7) !important;
        border: 1px solid rgba(0, 255, 204, 0.25) !important;
        border-radius: 8px !important;
        color: #FFFFFF !important;
        margin-bottom: 14px !important;
    }
    summary {
        color: #FFFFFF !important;
        font-weight: 600 !important;
        font-size: 0.95rem !important;
        padding: 8px 12px !important;
    }
</style>
""", unsafe_allow_html=True)


# ===========================================================================
# REQUIRED CUSTOM FUNCTIONS
# ===========================================================================

def retrieve_data(id):
    """
    Retrieve protein sequence data from NCBI Entrez database based on given ID.
    
    Parameters:
        id (str): Uniprot ID / Accession Number
        
    Returns:
        SeqRecord: Biopython SeqRecord object containing sequence and metadata
    """
    Entrez.email = "tanhuishan@graduate.utm.my"
    clean_id = id.strip()
    handle = Entrez.efetch(db='protein', id=clean_id, rettype='fasta', retmode='text')
    record = SeqIO.read(handle, 'fasta')
    handle.close()
    return record


def get_basic_analysis(sequence):
    """
    Analyze given protein sequence to calculate basic physicochemical properties.
    
    Parameters:
        sequence (str or Seq): Protein amino acid sequence
        
    Returns:
        tuple: (length, amino_acid_composition, molecular_weight, isoelectric_point)
    """
    clean_seq = str(sequence).strip().replace("\n", "").replace(" ", "").upper()
    seq_analysis = ProteinAnalysis(clean_seq)
    
    seq_length = len(clean_seq)
    aa_composition = seq_analysis.count_amino_acids()
    mol_weight = seq_analysis.molecular_weight()
    pI = seq_analysis.isoelectric_point()
    
    return seq_length, aa_composition, mol_weight, pI


# Amino Acid Lookup for detailed breakdown
AA_NAMES = {
    'A': 'Alanine (Ala)', 'C': 'Cysteine (Cys)', 'D': 'Aspartic Acid (Asp)',
    'E': 'Glutamic Acid (Glu)', 'F': 'Phenylalanine (Phe)', 'G': 'Glycine (Gly)',
    'H': 'Histidine (His)', 'I': 'Isoleucine (Ile)', 'K': 'Lysine (Lys)',
    'L': 'Leucine (Leu)', 'M': 'Methionine (Met)', 'N': 'Asparagine (Asn)',
    'P': 'Proline (Pro)', 'Q': 'Glutamine (Gln)', 'R': 'Arginine (Arg)',
    'S': 'Serine (Ser)', 'T': 'Threonine (Thr)', 'V': 'Valine (Val)',
    'W': 'Tryptophan (Trp)', 'Y': 'Tyrosine (Tyr)'
}


# ===========================================================================
# SIDEBAR NAVIGATION (Dark Biotech Minimalist)
# ===========================================================================
with st.sidebar:
    st.markdown('<h2 style="color: #FFFFFF; font-weight: 700; margin-top: 10px;">Bioinformatics II</h2>', unsafe_allow_html=True)
    st.markdown('<div style="color: #CBD5E1; font-size: 0.85rem; margin-bottom: 12px;">Lab Sequence Retrieval & Analysis Dashboard</div>', unsafe_allow_html=True)
    st.markdown("<hr style='border-color: rgba(0, 255, 204, 0.28);'>", unsafe_allow_html=True)
    
    selected_lab = st.radio("Lab", ["Lab 1"])


# ===========================================================================
# MAIN PAGE CONTENT
# ===========================================================================

if selected_lab == "Lab 1":
    # Header Section
    st.markdown('<div class="biotech-title">Lab 1 - Tan Hui Shan</div>', unsafe_allow_html=True)
    st.markdown('<div class="biotech-subtitle">Retrieval and Physicochemical Analysis of Protein Sequences via NCBI Entrez</div>', unsafe_allow_html=True)

    # Initialize empty session state so search bar starts BLANK on initial load
    if 'protein_id' not in st.session_state:
        st.session_state['protein_id'] = ""

    # Query Input Section (Rendered directly on page without enclosing rectangle)
    st.markdown('<div class="section-heading" style="margin-top: 0;">Query Parameters</div>', unsafe_allow_html=True)
    protein_id = st.text_input('Enter Uniprot ID', value=st.session_state['protein_id'], placeholder="e.g. P04637, P69905, P01308")
    retrieve = st.button('Retrieve & Analyze Sequence', type="primary")

    st.markdown("<hr style='border-color: rgba(0, 255, 204, 0.28); margin: 24px 0;'>", unsafe_allow_html=True)

    # Execution Flow: Runs when Retrieve button is clicked
    if retrieve:
        if protein_id.strip() != "":
            try:
                with st.spinner(f"Accessing NCBI Entrez repository for accession '{protein_id.strip()}'..."):
                    record = retrieve_data(protein_id.strip())
                    seq_len, aa_comp, mol_wt, pi_val = get_basic_analysis(record.seq)

                st.success("Sequence successfully retrieved and analyzed from database.")

                # -----------------------------------------------------------
                # SECTION 1: Protein Record Information (Glassmorphism card)
                # -----------------------------------------------------------
                st.markdown('<div class="section-heading">1. Protein Record Information</div>', unsafe_allow_html=True)
                
                st.markdown(f"""
                <div class="glass-card">
                    <div style="display: grid; grid-template-columns: 1fr 2fr; gap: 20px;">
                        <div>
                            <span class="field-label">Record Identifier</span>
                            <div class="monospace-text" style="font-size: 1.15rem; color: #00FFCC; font-weight: 700; margin-top: 4px; padding: 6px 10px;">{record.id}</div>
                            <span class="field-label" style="margin-top: 14px;">Record Name</span>
                            <div class="monospace-text" style="font-size: 1rem; color: #FFFFFF; margin-top: 4px; padding: 6px 10px;">{record.name}</div>
                        </div>
                        <div>
                            <span class="field-label">Description</span>
                            <div style="font-size: 1rem; color: #FFFFFF; line-height: 1.55; margin-top: 4px; background: rgba(5, 8, 17, 0.7); padding: 14px; border-radius: 6px; border: 1px solid rgba(0, 255, 204, 0.22); font-weight: 400;">{record.description}</div>
                        </div>
                    </div>
                </div>
                """, unsafe_allow_html=True)

                # -----------------------------------------------------------
                # SECTION 2: Genetic Sequence Data (Dropdown with full wrapped sequence)
                # -----------------------------------------------------------
                st.markdown('<div class="section-heading">2. Genetic Sequence Data</div>', unsafe_allow_html=True)
                
                # Format sequence in standard 60-character lines so all characters are directly visible
                seq_raw = str(record.seq)
                wrapped_sequence = "\n".join([seq_raw[i:i+60] for i in range(0, len(seq_raw), 60)])

                with st.expander("View Full Sequence & FASTA Record", expanded=False):
                    st.markdown(f'<div class="monospace-text" style="padding: 10px 14px; margin-bottom: 8px; color: #00FFCC; font-weight: 700;">>{record.description}</div>', unsafe_allow_html=True)
                    st.code(wrapped_sequence, language="text")

                with st.expander("View Raw Biopython Record Object", expanded=False):
                    st.code(str(record), language="text")

                # -----------------------------------------------------------
                # SECTION 3: Basic Protein Analysis (Combined Physicochemical Metrics + Composition)
                # -----------------------------------------------------------
                st.markdown('<div class="section-heading">3. Basic Protein Analysis</div>', unsafe_allow_html=True)
                
                # High-visibility clear stat cards for sequence length, molecular weight, and pI
                c1, c2, c3 = st.columns(3)
                with c1:
                    st.markdown(f"""
                    <div class="stat-card">
                        <div class="stat-label">Sequence Length</div>
                        <div class="stat-value">{seq_len} aa</div>
                        <div class="stat-sub">Total amino acid residues</div>
                    </div>
                    """, unsafe_allow_html=True)
                with c2:
                    st.markdown(f"""
                    <div class="stat-card">
                        <div class="stat-label">Molecular Weight</div>
                        <div class="stat-value">{mol_wt:,.2f} Da</div>
                        <div class="stat-sub">{(mol_wt / 1000):,.2f} kDa (Kilodaltons)</div>
                    </div>
                    """, unsafe_allow_html=True)
                with c3:
                    st.markdown(f"""
                    <div class="stat-card">
                        <div class="stat-label">Isoelectric Point</div>
                        <div class="stat-value">{pi_val:.3f} pI</div>
                        <div class="stat-sub">pH at net electrical neutrality</div>
                    </div>
                    """, unsafe_allow_html=True)

                st.markdown("<br>", unsafe_allow_html=True)

                # Amino Acid Composition Breakdown (Inside Basic Protein Analysis)
                df_composition = pd.DataFrame([
                    {
                        "Code": aa,
                        "Residue": AA_NAMES.get(aa, f"Residue {aa}"),
                        "Count": count,
                        "Percentage": round((count / seq_len) * 100, 2)
                    }
                    for aa, count in sorted(aa_comp.items())
                ])

                chart_col, table_col = st.columns([3, 2])

                # Rendered directly without enclosing rectangle cards
                with chart_col:
                    st.markdown('<div style="color: #00FFCC; font-size: 1.05rem; font-weight: 700; margin-bottom: 10px;">Residue Frequency Distribution</div>', unsafe_allow_html=True)
                    
                    # Create futuristic neon plot using Plotly Dark
                    fig = go.Figure()
                    fig.add_trace(go.Bar(
                        x=df_composition["Code"],
                        y=df_composition["Count"],
                        text=df_composition["Count"],
                        textposition="outside",
                        textfont=dict(color="#FFFFFF", size=11),
                        marker=dict(
                            color=df_composition["Count"],
                            colorscale=[[0, '#00A389'], [0.5, '#00FFCC'], [1.0, '#38BDF8']],
                            line=dict(color='#00FFCC', width=1.5)
                        ),
                        hovertemplate="<b>Residue:</b> %{x}<br><b>Frequency:</b> %{y}<extra></extra>"
                    ))
                    
                    fig.update_layout(
                        template="plotly_dark",
                        paper_bgcolor="rgba(0,0,0,0)",
                        plot_bgcolor="rgba(6, 10, 20, 0.7)",
                        xaxis=dict(
                            title=dict(text="Amino Acid Code", font=dict(color="#FFFFFF", size=13)),
                            gridcolor="rgba(0, 255, 204, 0.12)",
                            tickfont=dict(color="#FFFFFF", size=11)
                        ),
                        yaxis=dict(
                            title=dict(text="Abundance Count", font=dict(color="#FFFFFF", size=13)),
                            gridcolor="rgba(0, 255, 204, 0.12)",
                            tickfont=dict(color="#FFFFFF", size=11)
                        ),
                        margin=dict(l=45, r=20, t=25, b=45),
                        height=360
                    )
                    
                    st.plotly_chart(fig, use_container_width=True)

                with table_col:
                    st.markdown('<div style="color: #00FFCC; font-size: 1.05rem; font-weight: 700; margin-bottom: 10px;">Detailed Abundance Table</div>', unsafe_allow_html=True)
                    st.dataframe(
                        df_composition,
                        use_container_width=True,
                        hide_index=True,
                        height=360
                    )

            except Exception as e:
                st.error(f"Error accessing database or computing sequence analysis: {e}")
                st.info("Please verify the Uniprot ID and ensure internet connectivity to NCBI Entrez.")
        else:
            st.warning('Please enter a valid Uniprot ID')

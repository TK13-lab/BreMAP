"""Optional plotting of saved results; no fitting, relabeling or outcome selection."""
import numpy as np
import pandas as pd
import json


def plot_result(result, kind, *, x=None, y=None, annotations=None, group_by=None, save=None):
    """Return a Figure for programs, affinity, coverage or confusion.

    Programs show scored cells, including rejected cells; affinity is a cell-level or explicitly grouped mean
    heatmap, not cell abundance. Confusion uses candidate states for all joined
    cells, including rejected ones. Coverage always uses all query cells.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError('Install bremap[plot] to use visualization') from exc
    cells = result.table('mapping')
    if not cells.cell_id.is_unique:
        raise ValueError('Duplicate mapping cell IDs')
    fig, ax = plt.subplots(figsize=(7, 4.6), layout='constrained')
    if kind == 'coverage':
        groups = cells.groupby('donor_id', sort=True).accepted.agg(['sum','count'])
        coverage = groups['sum']/groups['count']
        ax.barh(groups.index.astype(str), coverage, color='#0072B2')
        for i, (_, row) in enumerate(groups.iterrows()):
            ax.text(.02, i, f"{int(row['sum'])}/{int(row['count'])}", va='center', color='white')
        ax.set(xlim=(0,1.02), xlabel='Accepted fraction (all query cells)', ylabel='Donor')
    elif kind == 'programs':
        programs = result.table('programs')
        if x not in programs or y not in programs or x == y or x in ('cell_id','donor_id') or y in ('cell_id','donor_id'):
            plt.close(fig)
            raise ValueError('Specify two different numeric program columns as x and y')
        data = cells.merge(programs, on=['cell_id','donor_id'], validate='one_to_one')
        if len(data) != len(cells):
            plt.close(fig)
            raise ValueError('Program/mapping cells do not match')
        colors = plt.get_cmap('tab10')
        for i, (state, group) in enumerate(data.groupby('candidate_state')):
            valid = group[[x,y]].notna().all(axis=1)
            ax.scatter(group.loc[valid,x], group.loc[valid,y], color=colors(i%10), label=state, s=40)
            rejected = valid & ~group.accepted
            ax.scatter(group.loc[rejected,x],group.loc[rejected,y], marker='x',s=120,color='black')
        ax.scatter([],[],marker='x',s=80,color='black',label='Rejected')
        units = json.loads((result.path/'workflow.json').read_text()).get('program_units','mean log1p(CP10k)')
        ax.set(xlabel=f'{x} program ({units})', ylabel=f'{y} program ({units})')
        ax.legend(loc='upper center',bbox_to_anchor=(.5,1.22),fontsize=9,ncol=2,frameon=False)
    elif kind == 'affinity':
        columns = [c for c in cells if c.startswith('conditional_affinity::')]
        if group_by is None:
            table = cells.sort_values(['candidate_state','cell_id']).set_index('cell_id')[columns]
        elif group_by in ('donor_id', 'candidate_state'):
            table = cells.groupby(group_by)[columns].mean()
        else:
            plt.close(fig)
            raise ValueError('group_by must be None, donor_id or candidate_state')
        fig.set_size_inches(7, min(14, max(4.6, len(table)*.24)))
        image = ax.imshow(table.to_numpy(),vmin=0,vmax=1,cmap='Blues',aspect='auto')
        ax.set_xticks(range(len(columns)),[c.split('::',1)[1] for c in columns],rotation=25,ha='right')
        ax.set_yticks(range(len(table)),table.index.astype(str),fontsize=8)
        ax.set(ylabel='Cell' if group_by is None else group_by, xlabel='Reference state')
        fig.colorbar(image,ax=ax,label='Conditional affinity' if group_by is None else 'Mean conditional affinity (scored cells)')
    elif kind == 'measurement':
        from .usability import compare_measurement
        if annotations is None:
            plt.close(fig)
            raise ValueError('Separate measurements required')
        joined, diagnostics = compare_measurement(cells, annotations, state=x, measurement=y)
        score = 'conditional_affinity::'+x
        valid = np.isfinite(joined[[score,y]]).all(axis=1)
        for accepted, label, color, marker in [(True,'Accepted','#0072B2','o'),(False,'Unassigned','#D55E00','x')]:
            block=joined.loc[valid & joined.accepted.eq(accepted)]
            ax.scatter(block[score],block[y],s=22,alpha=.65,color=color,marker=marker,label=f'{label} (n={len(block)})')
        ax.set(xlim=(-.02,1.02),xlabel=f'Reference affinity: {x}',ylabel=y)
        ax.legend(frameon=False)
    elif kind == 'confusion':
        if annotations is None:
            plt.close(fig)
            raise ValueError('Separate annotations required')
        labels = pd.read_csv(annotations) if not isinstance(annotations,pd.DataFrame) else annotations
        data = cells.merge(labels[['cell_id','donor_id','label']],on=['cell_id','donor_id'],validate='one_to_one')
        if len(data)!=len(cells) or data.label.isna().any():
            plt.close(fig)
            raise ValueError('Annotations must cover all mapped cells with matching donors')
        table = pd.crosstab(data.label,data.candidate_state)
        image=ax.imshow(table.to_numpy(),cmap='Blues',aspect='auto',vmin=0)
        for (i,j),value in np.ndenumerate(table.to_numpy()):
            ax.text(j,i,str(value),ha='center',va='center',color='white' if value>table.to_numpy().max()/2 else 'black')
        ax.set_xticks(range(len(table.columns)),table.columns,rotation=25,ha='right')
        ax.set_yticks(range(len(table)),table.index)
        ax.set(xlabel='Candidate state (all query cells)',ylabel='Separate annotation')
        fig.colorbar(image,ax=ax,label='Cells')
    else:
        plt.close(fig)
        raise ValueError('kind must be programs, affinity, coverage, measurement or confusion')
    if save is not None:
        fig.savefig(save,dpi=160,bbox_inches='tight')
    return fig


def plot_donor_correlations(data, *, baseline='rna_program_rho', reference='spearman_rho',
                            donor_column='donor_id',
                            labels=('Naive RNA program', 'BReMAP naive affinity'), save=None):
    """Plot paired within-donor correlations from a DataFrame or CSV path.

    Each row must identify one donor and two correlations with the same measured
    target on matched eligible cells. This function does not compute correlations
    or verify how their targets/folds were obtained. Values must be in [-1, 1];
    NaNs remain visible as missing, never converted to zero or silently dropped.
    Input row order is preserved, top to bottom. No significance test is performed.

    Parameters
    ----------
    data : pandas.DataFrame or path-like
        One row per donor. Defaults match the real-CD4 donor_validation.csv.
    baseline, reference : str
        Distinct numeric correlation columns, interpreted on the same scale.
    donor_column : str
        Unique, nonempty donor identifiers.
    labels : pair of str
        Baseline and reference legend labels, respectively.
    save : path-like, optional
        PNG/PDF/SVG destination; otherwise display the returned Figure in a notebook.

    Returns
    -------
    matplotlib.figure.Figure
        Customizable figure; plotting does not change input data.
    """
    table = data.copy() if isinstance(data, pd.DataFrame) else pd.read_csv(data, dtype={donor_column:str})
    if baseline == reference or donor_column in (baseline, reference):
        raise ValueError('Specify distinct donor and correlation columns')
    if table.empty or not set([donor_column, baseline, reference]) <= set(table):
        raise ValueError('Expected nonempty table with donor and both correlation columns')
    donors = table[donor_column]
    if donors.isna().any() or donors.astype(str).str.strip().eq('').any() or donors.astype(str).duplicated().any():
        raise ValueError('Donor identifiers must be unique and nonempty')
    if isinstance(labels, str) or len(labels) != 2 or any(not isinstance(v,str) or not v.strip() for v in labels):
        raise ValueError('Provide two nonempty legend labels')
    scores = table[[baseline, reference]].apply(pd.to_numeric, errors='raise').to_numpy(dtype=float, na_value=np.nan)
    if np.isinf(scores).any() or (np.abs(scores[~np.isnan(scores)]) > 1).any():
        raise ValueError('Correlations must be in [-1, 1]; use NaN for unavailable values')
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise ImportError('Install bremap[plot] to use visualization') from exc
    fig, ax = plt.subplots(figsize=(7, max(3.5, len(table)*.42+1.4)), layout='constrained')
    positions = np.arange(len(table))
    complete = ~np.isnan(scores).any(axis=1)
    for i in positions[complete]:
        ax.plot(scores[i], [i,i], color='#999999', zorder=1)
    for j, color in enumerate(['#777777', '#0072B2']):
        valid = ~np.isnan(scores[:,j])
        ax.scatter(scores[valid,j], positions[valid], label=labels[j], color=color, s=55, zorder=2)
    for i in positions[~complete]:
        missing = ', '.join(labels[j] for j in range(2) if np.isnan(scores[i,j]))
        ax.text(1.02, i, 'Missing: '+missing, transform=ax.get_yaxis_transform(),
                va='center', fontsize=9, color='#666666')
    ax.set_yticks(positions, donors.astype(str))
    ax.set(xlim=(-1,1), xlabel='Within-donor correlation with measured target', ylabel='Donor')
    ax.invert_yaxis()
    ax.axvline(0, color='#bbbbbb', lw=1, zorder=0)
    ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(.5,1.2))
    if save is not None:
        fig.savefig(save, dpi=160, bbox_inches='tight')
    return fig

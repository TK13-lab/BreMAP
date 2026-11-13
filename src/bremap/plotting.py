"""Optional plotting of saved results; no fitting, relabeling or outcome selection."""
import numpy as np
import pandas as pd


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
        ax.set(xlabel=f'{x} program (mean log1p CP10k)', ylabel=f'{y} program (mean log1p CP10k)')
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
        raise ValueError('kind must be programs, affinity, coverage or confusion')
    if save is not None:
        fig.savefig(save,dpi=160,bbox_inches='tight')
    return fig

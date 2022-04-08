def serialconnection():
    spacer_scaffold = pd.read_csv('spacer_scaffold_feature.csv')
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 33) & (
        spacer_scaffold.Pos_B <= 36), 'Structure'] = 'TL'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 54) & (
        spacer_scaffold.Pos_B <= 58), 'Structure'] = 'SL1'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 73) & (
        spacer_scaffold.Pos_B <= 76), 'Structure'] = 'SL2'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 88) & (
        spacer_scaffold.Pos_B <= 90), 'Structure'] = 'SL3'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 21) & (
        spacer_scaffold.Pos_B <= 32), 'Structure'] = 'R'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 37) & (
        spacer_scaffold.Pos_B <= 49), 'Structure'] = 'AR'
    spacer_scaffold.loc[(spacer_scaffold.Pos_B >= 63) & (
        spacer_scaffold.Pos_B <= 67), 'Structure'] = 'LR'
    spacer_scaffold.Structure.fillna('NS', inplace=True)
    spacer_scaffold.to_csv("Structural_annotation.csv", index=False)

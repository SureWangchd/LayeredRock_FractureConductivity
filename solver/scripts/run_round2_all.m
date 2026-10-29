function run_round2_all
% Reuses successful cases only if the complete saved parameter record matches.
groups={'smoke','screen','refine','mechanism','wave','contrast',...
        'structural_controls','structural_mechanism','tile_wave'};
for k=1:numel(groups),run_round2(groups{k});end
verify_round2_models;
enrich_round2_exchange;
enrich_round2_topology;
end

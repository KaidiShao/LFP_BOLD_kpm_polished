function dataset_id = get_dataset_id_from_stem(dataset_stem)
%GET_DATASET_ID_FROM_STEM Resolve the canonical dataset id from a file stem.

if nargin < 1 || isempty(dataset_stem)
    error('dataset_stem is required.');
end

dataset_stem = char(string(dataset_stem));
switch lower(dataset_stem)
    case 'e10fv1'
        dataset_id = 'E10.fV1';
    case 'e10gb1'
        dataset_id = 'E10.gb1';
    case 'e10gh1'
        dataset_id = 'E10.gH1';
    case 'e10gw1'
        dataset_id = 'E10.gW1';
    case 'e10aw1'
        dataset_id = 'E10.aW1';
    case 'e10bv1'
        dataset_id = 'E10.bv1';
    case 'e10ea1'
        dataset_id = 'E10.eA1';
    case 'f12m01'
        dataset_id = 'F12.m01';
    otherwise
        dataset_id = dataset_stem;
end
end

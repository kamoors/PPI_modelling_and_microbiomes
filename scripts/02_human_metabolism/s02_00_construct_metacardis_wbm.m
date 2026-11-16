function s02_00_construct_metacardis_wbm(hostModelDir, microbiomeModelDir, metadataFile, outputDir)
% SPDX-License-Identifier: PolyForm-Noncommercial-1.0.0
% Project-specific invocation of the external Persephone model-assembly pipeline.
% Install the study's Persephone/COBRA version separately; no pipeline code is bundled.
% Only input/output paths differ from the supervisor-provided source record.
% Host models, microbiome models and participant metadata are not distributed.

load(fullfile(hostModelDir, 'Harvey_1_03d.mat'));
modelH_male = male;

load(fullfile(hostModelDir, 'Harvetta_1_03d.mat'));
modelH_female = female;

%% Paths
MicrobiomeModels = microbiomeModelDir;
savePath = outputDir;
excelFile = metadataFile;

if ~exist(savePath, 'dir')
    mkdir(savePath);
end

%% Read metadata
T = readtable(excelFile);

files = dir(fullfile(MicrobiomeModels, '*.mat'));

for i = 1:length(files)
%for i = 1:5
    fileName = files(i).name;
    filePath = fullfile(MicrobiomeModels, fileName);

    % Extract ID
    tokens = regexp(fileName, 'microbiota_model_diet_(.+)\.mat', 'tokens');
    if isempty(tokens)
        fprintf('Skipping (pattern mismatch): %s\n', fileName);
        continue;
    end
    sampleID = tokens{1}{1};

    % Match metadata
    idx = strcmp(T.ID, sampleID);
    if ~any(idx)
        fprintf('No metadata for: %s\n', sampleID);
        continue;
    end

    gender = T.Gender{find(idx,1)};

    % Determine prefix + host model
    if strcmpi(gender, 'Male')
        modelH = modelH_male;
        prefix = 'Male_';
    elseif strcmpi(gender, 'Female')
        modelH = modelH_female;
        prefix = 'Female_';
    else
        fprintf('Unknown gender for %s\n', sampleID);
        continue;
    end

    % Expected output file
    saveName = fullfile(savePath, [prefix fileName]);

    % ---- SKIP if already processed ----
    if exist(saveName, 'file')
        fprintf('Skipping (already exists): %s\n', saveName);
        continue;
    end

    % Load microbiome model
    data = load(filePath);
    fn = fieldnames(data);
    modelM = data.(fn{1});

    % Combine
    try
        %modelHM = combineHarveyMicrotiota(modelH, modelM);
        modelHM = createMWBM(modelM,modelH, 'EUAverageDiet');
    catch ME
        fprintf('Error processing %s: %s\n', fileName, ME.message);
        continue;
    end

    %% Don't do diet adjustments

    % Save
    save(saveName, 'modelHM');

    fprintf('Processed: %s (%s)\n', fileName, gender);
end

end

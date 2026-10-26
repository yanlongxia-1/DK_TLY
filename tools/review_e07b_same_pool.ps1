param(
    [string]$CsvPath = (Join-Path $PSScriptRoot '..\archive\v0.2_fix\E07B_PhysicsFVIM_MVMD_v1\results_E07B\E07B_train_candidates.csv')
)

$ErrorActionPreference = 'Stop'
$resolved = (Resolve-Path -LiteralPath $CsvPath).Path
$rows = Import-Csv -LiteralPath $resolved | Where-Object { [int]$_.valid -eq 1 }
if (-not $rows) {
    throw "没有可用候选：$resolved"
}

$entropyBest = $rows | Sort-Object { [double]$_.entropyFitness } | Select-Object -First 1
$jointBest = $rows | Sort-Object { [double]$_.jointFitness } | Select-Object -First 1

[pscustomobject]@{
    Objective = 'entropyFitness'
    CandidateCount = $rows.Count
    K = [int]$entropyBest.K
    Alpha = [double]$entropyBest.alpha
    Fitness = [double]$entropyBest.entropyFitness
    OriginalRole = $entropyBest.role
}

[pscustomobject]@{
    Objective = 'jointFitness'
    CandidateCount = $rows.Count
    K = [int]$jointBest.K
    Alpha = [double]$jointBest.alpha
    Fitness = [double]$jointBest.jointFitness
    OriginalRole = $jointBest.role
}

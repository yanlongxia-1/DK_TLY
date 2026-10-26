$ErrorActionPreference = 'Stop'

$sources = @(
    'E:\dak_data\Desktop\实验\jyk_Test\DualGyro_FVIM_MVMD_Lab_v0.1',
    'E:\dak_data\Desktop\实验\jyk_Test\DualGyro_FVIM_MVMD_Lab_v0.2_fix'
)

$expected = @{
    'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-1-2\W3-10.TXT' = '590E5C79393948D50B4B83D14CE4A289CFA25D664A00E4C5B55CF6E6DBF425A3'
    'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-1-2\W3-20.TXT' = '9791C60F5462510D218C1CF236C74E4D0F4A4CFD0969EA9B8683B1A8D7BFD41A'
    'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-3\Z3-1.TXT' = '2BD88AA312EE2662550FFACE6EBC4DF454C53C5D9C4BA35B7DF728E22FD33F92'
    'E:\dak_data\Desktop\实验\jyk_Test\实验室课题组数据实验\20200824_data\20200824陀螺试验\20200824-3\Z3-2.TXT' = 'F18D265D32B6C4581A63EAE8F7148966C3F96F990314E58B79C268FC72278833'
}

foreach ($source in $sources) {
    if (-not (Test-Path -LiteralPath $source)) {
        throw "缺少只读来源目录：$source"
    }
    Write-Host "来源目录存在：$source"
}

foreach ($entry in $expected.GetEnumerator()) {
    if (-not (Test-Path -LiteralPath $entry.Key)) {
        throw "缺少关键真实数据：$($entry.Key)"
    }
    $actual = (Get-FileHash -LiteralPath $entry.Key -Algorithm SHA256).Hash
    if ($actual -ne $entry.Value) {
        throw "SHA-256 不匹配：$($entry.Key)`n期望：$($entry.Value)`n实际：$actual"
    }
    Write-Host "哈希通过：$($entry.Key)"
}

Write-Host '来源与四个关键真实数据文件核验通过。'

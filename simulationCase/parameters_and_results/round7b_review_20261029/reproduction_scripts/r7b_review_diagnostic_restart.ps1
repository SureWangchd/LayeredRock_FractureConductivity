$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
$out=Join-Path $root 'runs/round7b_review_20261029'
$state=Join-Path $out 'tolerance_checks/S1_M3_water_0p18_T_background_tol_1em11'
if(Test-Path -LiteralPath "$state/comparison.json"){throw 'First state completed; do not interrupt it.'}
$request=Get-Item -LiteralPath "$state/request_before_solve.json"
$elapsed=((Get-Date)-$request.LastWriteTime).TotalSeconds
if($elapsed -lt 1500){throw "25-minute diagnostic checkpoint has not elapsed: $elapsed seconds."}
$ownedId=[int](Get-Content -LiteralPath "$out/server.pid")
$owned=Get-CimInstance Win32_Process -Filter "ProcessId=$ownedId"
if(-not $owned -or $owned.Name -ne 'comsolmphserver.exe' -or $owned.CommandLine -notmatch '-port 2038'){throw 'Review server ownership mismatch.'}
$attempt=Join-Path $out 'attempts/zero_start_20261029_221636'
New-Item -ItemType Directory -Path $attempt -Force | Out-Null
[ordered]@{status='STOPPED_TO_ENABLE_DIAGNOSTICS_AND_SAVED_INITIAL_GUESS';time_local=(Get-Date -Format s);elapsed_since_request_seconds=$elapsed;server_pid=$ownedId;reason='No return after the 25-minute checkpoint and no per-linear-iteration log. Repeat the same fixed mesh, materials, frozen kappa and 1e-11 tolerance with a verified saved initial field and detailed logs.';interpretation='This active termination does not establish nonconvergence or roundoff stagnation.';source_outputs_overwritten=$false} | ConvertTo-Json | Set-Content -LiteralPath "$attempt/ACTIVE_TERMINATION.json"
Stop-Process -Id $ownedId
Start-Sleep -Seconds 6
Copy-Item -LiteralPath $state -Destination "$attempt/T_background_attempt" -Recurse
Copy-Item -LiteralPath "$out/tolerance_checks/S1_M3_water_0p18_pair_four_states" -Destination "$attempt/pair_attempt" -Recurse
Get-ChildItem -LiteralPath $out -File | Where-Object {$_.Extension -eq '.log' -or $_.Name -like 'resources_*.csv' -or $_.Name -like 'FAILURE_*.json'} | Copy-Item -Destination $attempt
Write-Output "Archived first attempt: $attempt"
& "$PSScriptRoot/r7b_review_launch.ps1"

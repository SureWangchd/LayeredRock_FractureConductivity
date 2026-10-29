param([ValidateRange(1,32)][int]$Threads=26)
$ErrorActionPreference='Stop'
$projectRoot=Split-Path -Parent $PSScriptRoot
$out=Join-Path $projectRoot 'runs/round7b_review_20261029'
New-Item -ItemType Directory -Path $out -Force | Out-Null
$listener=Get-NetTCPConnection -LocalPort 2038 -State Listen -ErrorAction SilentlyContinue
if($listener){throw 'Review port 2038 already in use; do not attach another job.'}
$stamp=Get-Date -Format yyyyMMdd_HHmmss
$exe='D:/Program Files/COMSOL/COMSOL63/Multiphysics/bin/win64/comsolmphserver.exe'
$server=Start-Process -FilePath $exe -ArgumentList "-np $Threads -port 2038 -login auto -multi on -silent -autosave off" -WindowStyle Hidden -RedirectStandardOutput "$out/server_$stamp.stdout.log" -RedirectStandardError "$out/server_$stamp.stderr.log" -PassThru
$server.Id | Set-Content -LiteralPath "$out/server.pid"
[ordered]@{time_local=(Get-Date -Format s);threads=$Threads;server_pid=$server.Id;port=2038;source_mesh_and_kappa_fixed=$true;wall_time_kill=$false} | ConvertTo-Json | Set-Content -LiteralPath "$out/launch_$stamp.json"
$ready=$false
for($attempt=0;$attempt -lt 50;$attempt++){
 $check=Get-NetTCPConnection -LocalPort 2038 -State Listen -ErrorAction SilentlyContinue
 if($check -and ($check.OwningProcess -contains $server.Id)){$ready=$true;break}
 if(-not (Get-Process -Id $server.Id -ErrorAction SilentlyContinue)){throw 'Review server exited.'}
 Start-Sleep -Seconds 1
}
if(-not $ready){throw 'Review server not ready.'}
$watch=Start-Process -FilePath powershell.exe -ArgumentList @('-NoProfile','-ExecutionPolicy','Bypass','-File',('"'+(Join-Path $PSScriptRoot 'r7b_review_watch.ps1')+'"')) -WindowStyle Hidden -PassThru
$watch.Id | Set-Content -LiteralPath "$out/resource_watch.pid"
$codePath=$PSScriptRoot.Replace([char]92,[char]47).Replace("'","''")
$safeRoot=$projectRoot.Replace([char]92,[char]47).Replace("'","''")
& 'E:/matlab/bin/matlab.exe' -batch "addpath('$codePath');r7b_review_run('$safeRoot')" *> "$out/launch_$stamp.log"
if($LASTEXITCODE -ne 0){throw "Review MATLAB exited $LASTEXITCODE; check preserved logs."}

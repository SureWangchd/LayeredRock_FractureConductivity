$ErrorActionPreference='Stop'
$root=Split-Path -Parent $PSScriptRoot
$out=Join-Path $root 'runs/round7b_review_20261029'
$ownedId=[int](Get-Content -LiteralPath "$out/server.pid")
$owned=Get-CimInstance Win32_Process -Filter "ProcessId=$ownedId"
if(-not $owned -or $owned.Name -ne 'comsolmphserver.exe' -or $owned.CommandLine -notmatch '-port 2038'){throw 'Review server ownership mismatch.'}
$log=Join-Path $out ('resources_'+(Get-Date -Format yyyyMMdd_HHmmss)+'.csv')
'time_local,free_memory_MB,server_working_set_MB,server_private_MB,cpu_seconds' | Set-Content -LiteralPath $log
$low=0
while(Get-Process -Id $ownedId -ErrorAction SilentlyContinue){
 $os=Get-CimInstance Win32_OperatingSystem
 $process=Get-Process -Id $ownedId -ErrorAction SilentlyContinue
 if(-not $process){break}
 $free=[double]$os.FreePhysicalMemory/1024
 '{0},{1:F1},{2:F1},{3:F1},{4:F1}' -f (Get-Date -Format 'yyyy-MM-dd HH:mm:ss'),$free,($process.WorkingSet64/1MB),($process.PrivateMemorySize64/1MB),$process.CPU | Add-Content -LiteralPath $log
 if($free -lt 5000){$low++}else{$low=0}
 if($low -ge 2){
  [ordered]@{status='stopped_for_memory_reserve';time_local=(Get-Date -Format s);free_MB=$free;server_pid=$ownedId} | ConvertTo-Json | Set-Content -LiteralPath "$out/resource_stop_$(Get-Date -Format yyyyMMdd_HHmmss).json"
  Stop-Process -Id $ownedId
  break
 }
 Start-Sleep -Seconds 5
}

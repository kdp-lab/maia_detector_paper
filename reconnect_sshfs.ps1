# Define SSHFS mount details
$driveLetter = "X:"
$remotePath = "\\sshfs.r\lrozanov@kdplab01.uchicago.edu"
$userName = "lrozanov"
$server = "kdplab01"
$password = "League_Of_Legends" # Consider using a more secure method for handling passwords

# Check if the drive is already mounted
if (!(Test-Path $driveLetter)) {
    Write-Host "Mount not found. Are you sure you're using the VPN? Attempting to reconnect..."

    # Attempt to mount the SSHFS drive
    $credential = "$userName:$password"
    $connectResult = net use $driveLetter $remotePath $credential /user:$userName

    if ($connectResult -like "*The command completed successfully.*") {
        Write-Host "Successfully reconnected the SSHFS mount."
    } else {
        Write-Host "Failed to reconnect the SSHFS mount."
    }
} else {
    Write-Host "SSHFS mount is currently available."
}

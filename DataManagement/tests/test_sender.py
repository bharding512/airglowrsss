import os
from unittest import mock
from unittest.mock import Mock, patch
from botocore.exceptions import ClientError

from airglow.Sender import send


def filter_system_commands(calls: list[list[str]], command:str) -> list[str]:
    # Filter sublists where the first element starts with "chmod"
    filtered_list = [sublist[0][0] for sublist in calls if
     sublist and isinstance(sublist[0][0], str) and sublist[0][0].startswith(command)]
    return filtered_list

@patch("airglow.Sender.os.system")
@patch("airglow.Sender.os.remove")
def test_send(mock_remove, mock_system):
    mock_system.return_value = 0
    send({"./sample_files"}, "airglow.illinois.edu:/home/airglow/Sending/", 10)
    chmod_calls = filter_system_commands(mock_system.call_args_list, "chmod")
    assert len(chmod_calls) == 2
    assert chmod_calls == ['chmod 774 ./sample_files/bar.tar.gzip', 'chmod 774 ./sample_files/foo.txt']

    scp_calls = filter_system_commands(mock_system.call_args_list, "scp")
    assert len(scp_calls) == 2
    assert scp_calls == ['scp ./sample_files/bar.tar.gzip airglow.illinois.edu:/home/airglow/Sending/', 'scp ./sample_files/foo.txt airglow.illinois.edu:/home/airglow/Sending/']

    assert mock_remove.call_count == 2
    removes = [call[0][0] for call in mock_remove.call_args_list]
    assert removes == [
        "./sample_files/bar.tar.gzip",
        "./sample_files/foo.txt"
    ]

@patch("airglow.Sender.os.system")
@patch("airglow.Sender.os.remove")
def test_send_ssh_err(mock_remove, mock_system):
    mock_system.side_effect = [0, -1, 0, 0]
    send({"./sample_files"}, "airglow.illinois.edu:/home/airglow/Sending/", 10)
    chmod_calls = filter_system_commands(mock_system.call_args_list, "chmod")
    assert len(chmod_calls) == 2

    scp_calls = filter_system_commands(mock_system.call_args_list, "scp")
    assert len(scp_calls) == 2

    assert mock_remove.call_count == 1
    removes = [call[0][0] for call in mock_remove.call_args_list]
    assert removes == [
        "./sample_files/foo.txt"
    ]

@patch("airglow.Sender.boto3.client")
@patch("airglow.Sender.os.system")
@patch("airglow.Sender.os.remove")
@mock.patch.dict(os.environ, {
    "AWS_ACCESS_KEY_ID": "test-key",
    "AWS_SECRET_ACCESS_KEY": "test-secret",
    "AWS_S3_ENDPOINT_URL": "http://localhost:4566",
    "DEST_BUCKET": "test-bucket"
})
def test_s3_upload(mock_remove, mock_system, mock_boto3):
    mock_boto3.return_value = Mock()
    mock_boto3.return_value.upload_file = Mock(
        side_effect=[None, ClientError]
    )

    mock_system.return_value = 0

    send({"./sample_files"}, "airglow.illinois.edu:/home/airglow/Sending/", 10)
    assert mock_boto3.return_value.upload_file.call_count == 2
    assert mock_boto3.return_value.upload_file.call_args_list[0].args == (
        "./sample_files/bar.tar.gzip",
        "test-bucket",
        "airglow/raw/./sample_files/bar.tar.gzip"
    )

    assert mock_system.call_count == 2 * 2  # 2 files, 2 commands each
    assert mock_remove.call_count == 1  # Only one file was uploaded successfully
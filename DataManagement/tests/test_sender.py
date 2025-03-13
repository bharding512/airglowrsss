from unittest.mock import patch

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
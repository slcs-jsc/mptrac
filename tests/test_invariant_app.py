import pytest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "projects", "web_runner"))

from app import app


@pytest.fixture
def client():
    app.config["TESTING"] = True
    with app.test_client() as client:
        yield client


@pytest.mark.parametrize("payload", [
    {"cmd": "echo hello; rm -rf /", "expected_absent": "rm"},
    {"cmd": "test$(whoami)", "expected_absent": None},
    {"cmd": "ls | cat /etc/passwd", "expected_absent": "root:"},
    {"cmd": "valid_command", "expected_absent": None},
])
def test_no_shell_injection_via_run_endpoint(payload, client):
    """Invariant: Shell metacharacters in user input must not result in arbitrary command execution."""
    response = client.post("/run", data={"cmd": payload["cmd"]})
    
    # The response should never contain evidence of injected command execution
    if payload["expected_absent"]:
        response_text = response.get_data(as_text=True)
        assert payload["expected_absent"] not in response_text, (
            f"Shell injection detected: found '{payload['expected_absent']}' in response "
            f"when submitting payload '{payload['cmd']}'"
        )
    
    # Regardless of input, server should not crash (500 is acceptable for bad input,
    # but we should not see signs of successful injection)
    assert response.status_code in (200, 400, 403, 404, 422, 500)
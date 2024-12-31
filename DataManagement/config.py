# config.py
import os
from dotenv import load_dotenv
from pathlib import Path

class Config:
    def __init__(self):
        # Find .env file from current directory up to root
        env_path = Path('.').absolute()
        while env_path != env_path.parent:
            if (env_path / '.env').exists():
                break
            env_path = env_path.parent
        
        # Load environment variables from .env file
        load_dotenv(env_path / '.env')
        
        # Get configuration values with defaults
        self.github_token = os.getenv('GITHUB_TOKEN')
        self.github_repo = os.getenv('GITHUB_REPO')
        
        # Validate required settings
        if not self.github_token:
            raise ValueError("GitHub token not found in environment variables")
        if not self.github_repo:
            raise ValueError("GitHub repository not found in environment variables")

# Create a singleton instance
config = Config()
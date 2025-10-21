"""
Simple script to convert ChatGPT conversation to Jupyter notebook cells
Run this in a Jupyter notebook cell
"""

import json
import re

def convert_chatgpt_to_notebook(conversation_text, notebook_title="ChatGPT Conversation"):
    """
    Convert ChatGPT conversation text to Jupyter notebook format
    
    Args:
        conversation_text (str): The full ChatGPT conversation text
        notebook_title (str): Title for the notebook
    
    Returns:
        dict: Jupyter notebook structure
    """
    
    cells = []
    
    # Pattern to match code blocks with optional language specification
    code_pattern = r'```(?:python|py|javascript|js|bash|sh|sql|r|matlab|markdown|md)?\n(.*?)```'
    
    # Split text by code blocks
    parts = re.split(code_pattern, conversation_text, flags=re.DOTALL)
    
    for i, part in enumerate(parts):
        part = part.strip()
        if not part or part.isspace():
            continue
            
        if i % 2 == 0:  # Text parts (even indices)
            # Add as markdown cell
            cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [part]
            })
        else:  # Code parts (odd indices)
            # Add as code cell
            cells.append({
                "cell_type": "code",
                "execution_count": None,
                "metadata": {},
                "outputs": [],
                "source": [part]
            })
    
    # Create notebook structure
    notebook = {
        "cells": cells,
        "metadata": {
            "kernelspec": {
                "display_name": "Python 3",
                "language": "python",
                "name": "python3"
            },
            "language_info": {
                "codemirror_mode": {
                    "name": "ipython",
                    "version": 3
                },
                "file_extension": ".py",
                "mimetype": "text/x-python",
                "name": "python",
                "nbconvert_exporter": "python",
                "pygments_lexer": "ipython3",
                "version": "3.8.0"
            }
        },
        "nbformat": 4,
        "nbformat_minor": 4
    }
    
    return notebook

# Example usage:
def save_notebook(notebook, filename):
    """Save notebook to file"""
    with open(filename, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)
    print(f"Notebook saved as: {filename}")

# Example:
# 1. Copy your ChatGPT conversation text
# 2. Paste it as a string variable
# 3. Run the conversion

# conversation = """
# Your ChatGPT conversation text here...
# """

# notebook = convert_chatgpt_to_notebook(conversation, "My ChatGPT Conversation")
# save_notebook(notebook, "my_conversation.ipynb")

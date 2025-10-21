#!/usr/bin/env python3
"""
Convert ChatGPT conversation to Jupyter notebook (.ipynb) format
"""

import json
import re
import argparse
from pathlib import Path

def parse_chatgpt_conversation(text):
    """
    Parse ChatGPT conversation text and extract code blocks and text.
    Returns a list of cells for the notebook.
    """
    cells = []
    
    # Split by common ChatGPT delimiters
    # Look for code blocks (```python, ```, etc.)
    code_pattern = r'```(?:python|py|javascript|js|bash|sh|sql|r|matlab)?\n(.*?)```'
    text_pattern = r'```(?:python|py|javascript|js|bash|sh|sql|r|matlab)?\n.*?```'
    
    # Find all code blocks
    code_blocks = re.findall(code_pattern, text, re.DOTALL)
    text_blocks = re.split(text_pattern, text, flags=re.DOTALL)
    
    # Process text and code blocks
    for i, text_block in enumerate(text_blocks):
        text_block = text_block.strip()
        if text_block and not text_block.isspace():
            # Add as markdown cell
            cells.append({
                "cell_type": "markdown",
                "metadata": {},
                "source": [text_block]
            })
        
        # Add corresponding code block if it exists
        if i < len(code_blocks):
            code_content = code_blocks[i].strip()
            if code_content:
                cells.append({
                    "cell_type": "code",
                    "execution_count": None,
                    "metadata": {},
                    "outputs": [],
                    "source": [code_content]
                })
    
    return cells

def create_notebook(cells, title="ChatGPT Conversation"):
    """
    Create a Jupyter notebook structure from cells.
    """
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

def main():
    parser = argparse.ArgumentParser(description='Convert ChatGPT conversation to Jupyter notebook')
    parser.add_argument('input_file', help='Input text file containing ChatGPT conversation')
    parser.add_argument('-o', '--output', help='Output notebook file (default: input_file.ipynb)')
    parser.add_argument('-t', '--title', default='ChatGPT Conversation', help='Notebook title')
    
    args = parser.parse_args()
    
    # Read input file
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file '{input_file}' not found")
        return
    
    with open(input_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Parse conversation
    cells = parse_chatgpt_conversation(content)
    
    if not cells:
        print("No content found to convert")
        return
    
    # Create notebook
    notebook = create_notebook(cells, args.title)
    
    # Determine output file
    if args.output:
        output_path = Path(args.output)
    else:
        output_path = input_path.with_suffix('.ipynb')
    
    # Write notebook
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=2, ensure_ascii=False)
    
    print(f"Notebook created: {output_path}")
    print(f"Total cells: {len(cells)}")

if __name__ == "__main__":
    main()








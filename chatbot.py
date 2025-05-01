#!/usr/bin/env python3
from openai import OpenAI
import sys

# Model aliases
MODEL_MAP = {
    "gpt3": "openai/gpt-3.5-turbo",
    "gpt4": "openai/gpt-4o-2024-11-20",
    "o3mh": "openai/o3-mini-high",
    "deepseek": "deepseek/deepseek-r1",
    "claude": "anthropic/claude-3.5-sonnet-20240620",
    "gemini": "google/gemini-2.0-flash-001"
}

def ask_chatbot(user_question, model_choice = 'gpt3'):
    model_choice = model_choice.lower()
    if model_choice not in MODEL_MAP:
        raise ValueError(f"Invalid model choice: {model_choice}")

    model = MODEL_MAP[model_choice]
    client = OpenAI(
        base_url="https://openrouter.ai/api/v1",
        api_key="YOUR_API_KEY"
    )

    completion = client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": user_question}]
    )

    return completion.choices[0].message.content

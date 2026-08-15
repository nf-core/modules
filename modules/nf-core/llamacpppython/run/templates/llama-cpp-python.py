#!/usr/bin/env python3

import argparse
import json
import os
import shlex
import sys

import llama_cpp


# Helper to create messages from a text file
def create_messages_from_textfile(textfile, system_prompt):
    try:
        with open(textfile, encoding="utf-8") as f:
            content = f.read()
        return [
            {"role": "system", "content": system_prompt.strip()},
            {"role": "user", "content": content.strip()},
        ]
    except Exception as e:
        print(f"Error reading text file '{textfile}': {e}", file=sys.stderr)
        sys.exit(1)


# Helper to load messages from JSON or fallback to text
def load_messages(messages_file, system_prompt):
    if not os.path.exists(messages_file):
        print(f"Messages file '{messages_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    try:
        with open(messages_file, encoding="utf-8") as f:
            content = f.read()
            try:
                return json.loads(content)
            except json.JSONDecodeError:
                return create_messages_from_textfile(messages_file, system_prompt)
    except Exception as e:
        print(f"Error opening messages file '{messages_file}': {e}", file=sys.stderr)
        sys.exit(1)


def llamacpp_python(
    messages_file,
    model_file,
    temperature=0.9,
    output="output.txt",
    verbose=False,
    context_size=2048,
    chat_format="chatml",
    seed=None,
):
    if not os.path.exists(model_file):
        print(f"Model file '{model_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    # Default system prompt
    system_prompt = "A chat between a curious user and an artificial intelligence assistant. The assistant gives helpful, detailed, and polite answers to the user's questions"

    messages_json = load_messages(messages_file, system_prompt)

    try:
        llm = llama_cpp.Llama(
            model_path=model_file,
            chat_format=chat_format,
            n_ctx=context_size,
            seed=seed,
        )
        response = llm.create_chat_completion(
            messages=messages_json,
            response_format={"type": "json_object"},
            temperature=temperature,
        )
    except Exception as e:
        print(f"Error running llama_cpp: {e}", file=sys.stderr)
        sys.exit(1)

    if not verbose:
        try:
            reply = response["choices"][0]["message"]["content"]
        except (KeyError, IndexError, TypeError):
            reply = response
        # Try to parse reply as JSON if it's a string
        if isinstance(reply, str):
            try:
                reply_json = json.loads(reply)
                if isinstance(reply_json, dict) and len(reply_json) == 1:
                    reply = next(iter(reply_json.values()))
                else:
                    reply = reply_json
            except Exception:
                pass  # Leave reply as string if not valid JSON
        elif isinstance(reply, dict) and len(reply) == 1:
            reply = next(iter(reply.values()))
    else:
        reply = response

    try:
        with open(output, "w", encoding="utf-8") as f:
            if isinstance(reply, str):
                f.write(reply)
            else:
                f.write(json.dumps(reply, indent=2))
        if verbose:
            print(f"Output written to {output}")
    except Exception as e:
        print(f"Error writing output file '{output}': {e}", file=sys.stderr)
        sys.exit(1)


def main(args_string=None):
    parser = argparse.ArgumentParser(description="Submit a process with model.")
    parser.add_argument("-s", "--messages", required=True, help="JSON message")
    parser.add_argument("-m", "--model", required=True, help="Model used")
    parser.add_argument("-t", "--temperature", default=0.9, type=float, help="Temperature")
    parser.add_argument("-o", "--output", default="output.txt", help="Output text")
    parser.add_argument("-c", "--context", default=2048, type=int, help="Context size")
    parser.add_argument("--chat_format", default="chatml", help="Chat format")
    parser.add_argument("--seed", default=None, type=int, help="Defined seed")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")

    args = parser.parse_args(shlex.split(args_string) if args_string is not None else None)
    llamacpp_python(
        messages_file=args.messages,
        model_file=args.model,
        temperature=args.temperature,
        output=args.output,
        verbose=args.verbose,
        context_size=args.context,
        chat_format=args.chat_format,
        seed=args.seed,
    )


def get_cuda_version():
    """Return CUDA runtime version (major.minor) or 'no CUDA available'."""
    import ctypes

    for major in (13, 12, 11):
        try:
            cudart = ctypes.CDLL(f"libcudart.so.{major}")
        except OSError:
            continue
        v = ctypes.c_int()
        if cudart.cudaRuntimeGetVersion(ctypes.byref(v)) == 0:
            return f"{v.value // 1000}.{(v.value % 1000) // 10}"
    return "no CUDA available"


def write_versions():
    versions = {
        "${task.process}": {
            "llama-cpp-python": llama_cpp.__version__,
            "cuda": get_cuda_version(),
        }
    }
    with open("versions.yml", "w", encoding="utf-8") as f:
        for process, pkgs in versions.items():
            f.write(f'"{process}":\\n')
            for pkg, ver in pkgs.items():
                f.write(f"    {pkg}: {ver}\\n")


if __name__ == "__main__":
    main("--model ${gguf_model} --messages ${prompt_file} --output ${prefix}.txt ${args}")
    write_versions()

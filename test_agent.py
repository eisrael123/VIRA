import asyncio
from claude_agent_sdk import query, ClaudeAgentOptions, tool
from claude_agent_sdk.types import AssistantMessage, ResultMessage, TextBlock

# 1. Define a tool that Claude can call to print consumption metrics.
# The SDK requires the tool decorator with clear name, description, and input schema.
@tool(
    name="report_token_usage_and_cost",
    description="Use this tool at the very end of your response to calculate and print estimated token usage and USD costs.",
    input_schema={
        "type": "object",
        "properties": {
            "input_tokens": {"type": "integer", "description": "Approximate number of input tokens used."},
            "output_tokens": {"type": "integer", "description": "Approximate number of output tokens generated."}
        },
        "required": ["input_tokens", "output_tokens"]
    }
)
async def report_token_usage_and_cost(input_tokens: int, output_tokens: int):
    # This executes locally on your machine when Claude calls it
    print(f"\n\n[Tool Executed] Reporting footprint: ~{input_tokens} Input / ~{output_tokens} Output tokens.")
    return {"status": "Logged successfully"}


async def run_test():
    prompt = "Write a short, single-paragraph poem about a cat running outside. At the very bottom, invoke the report_token_usage_and_cost tool."
    
    # Configure options. permission_mode="acceptEdits" allows autonomous execution
    # of tools without blocking for user terminal confirmations (y/n).
    options = ClaudeAgentOptions(
        allowed_tools=["report_token_usage_and_cost"],
        permission_mode="acceptEdits",
        system_prompt="You are a helpful test assistant. You must always use your custom tools when requested."
    )
    
    print("Sending prompt to Claude Agent SDK...")
    print("-" * 50)
    
    # query returns an async iterator streaming events across multiple turns
    async for message in query(prompt=prompt, options=options):
        
        # Capture streaming text chunks as Claude thinks/responds
        if isinstance(message, AssistantMessage):
            for block in message.content:
                if isinstance(block, TextBlock):
                    print(block.text, end="", flush=True)
                    
        # Capture the final loop summary containing exact client-side cost metrics
        elif isinstance(message, ResultMessage):
            print("\n" + "-" * 50)
            print("--- FINAL CLIENT-SIDE ACCOUNTING ---")
            print(f"Total Estimated Cost: ${message.total_cost_usd:.5f} USD")
            
if __name__ == "__main__":
    # Execute the async loop inside the terminal
    asyncio.run(run_test())
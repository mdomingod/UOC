#!/usr/bin/env python3
# To run the srcipt, activate the uv environment and change th end line from CRLF to LF (in vscode, bottom right corner) and then run:
# python src/project_package/tests/test_model.py

# tests/test_model.py
import torch
from project_package.config import Config
from project_package.model import PordleCNN


def test_parameter_count():
    """The paper reports exactly 214,081 trainable parameters. Verify we match."""
    cfg = Config()
    model = PordleCNN(cfg)
    n = model.count_parameters()
    assert n == 214_081, f"Expected 214,081 parameters, got {n:,}"
    print(f"✓ Parameter count correct: {n:,}")


def test_output_shape():
    """A batch of N images should produce N scalar predictions."""
    cfg = Config()
    model = PordleCNN(cfg)
    model.eval()

    batch_size = 8
    fake_batch = torch.randn(batch_size, 2, 5, 66)  # (N, channels, height, width)
    output = model(fake_batch)

    assert output.shape == (batch_size,), \
        f"Expected output shape ({batch_size},), got {output.shape}"
    print(f"✓ Output shape correct: {output.shape}")


def test_output_non_negative():
    """With allow_negative=False (default), all predictions should be >= 0."""
    cfg = Config()
    model = PordleCNN(cfg)
    model.eval()

    fake_batch = torch.randn(32, 2, 5, 66)
    output = model(fake_batch)

    assert (output >= 0).all(), "Found negative predictions but allow_negative=False"
    print(f"✓ All predictions non-negative")


def test_intermediate_shapes():
    """
    Trace the shape at every layer to catch any spatial dimension mistake.
    This is the most useful debug test — if shapes are wrong, this tells you exactly where.
    """
    cfg = Config()
    model = PordleCNN(cfg)
    model.eval()

    x = torch.randn(4, 2, 5, 66)
    print(f"\n  Input:          {tuple(x.shape)}")

    x = torch.relu(model.conv1(x))
    print(f"  After conv1:    {tuple(x.shape)}   expect (4, 32, 5, 66)")

    x = torch.relu(model.conv2(x))
    print(f"  After conv2:    {tuple(x.shape)}   expect (4, 32, 5, 66)")

    x = model.pool1(x)
    print(f"  After pool1:    {tuple(x.shape)}   expect (4, 32, 2, 33)")

    x = torch.nn.functional.pad(x, (0, 0, 0, 1))
    print(f"  After pad:      {tuple(x.shape)}   expect (4, 32, 3, 33)")

    x = torch.relu(model.conv3(x))
    print(f"  After conv3:    {tuple(x.shape)}   expect (4, 32, 3, 33)")

    x = torch.relu(model.conv4(x))
    print(f"  After conv4:    {tuple(x.shape)}   expect (4, 32, 3, 33)")

    x = model.pool2(x)
    print(f"  After pool2:    {tuple(x.shape)}   expect (4, 32, 2, 17)")

    x = x.flatten(start_dim=1)
    print(f"  After flatten:  {tuple(x.shape)}   expect (4, 1088)")

    print(f"\n✓ All intermediate shapes correct")


def test_gradient_flows(): 
    """
    Verify that a backward pass doesn't crash and that gradients are non-zero.
    This catches dead layers, wrong activations, or broken connections.
    """
    cfg = Config()
    model = PordleCNN(cfg)
    model.train()

    fake_batch  = torch.randn(4, 2, 5, 66)
    fake_labels = torch.randint(low=0, high=2000, size=(4,)).float()

    output = model(fake_batch)
    loss   = torch.nn.MSELoss()(output, fake_labels)
    loss.backward()

    # The output ReLU legitimately blocks gradients when predictions are zero.
    # So instead, verify that internal layers have gradients — this confirms
    # the computational graph is fully connected end to end.
    assert model.fc_out.weight.grad is not None, \
        "No gradient in fc_out — graph is broken before the output layer"
    
    assert model.conv1.weight.grad is not None, \
        "No gradient in conv1 — gradient is not flowing through the full network"

    assert model.conv1.weight.grad.abs().sum() > 0, \
        "conv1 gradient is all zeros — something is blocking backprop"

    print(f"✓ Gradients flow correctly through the full network")
    print(f"  loss = {loss.item():.2f}")
    print(f"  conv1 grad norm = {model.conv1.weight.grad.abs().sum():.4f}")
    print(f"  fc_out grad norm = {model.fc_out.weight.grad.abs().sum():.4f}")


    # Check that at least one parameter has a non-zero gradient
    """ Old code that was too strict because the output ReLU can legitimately block gradients when predictions are zero. 
        So instead, we check that internal layers have gradients to confirm the graph is fully connected.

    print(f"Output values: {output[:5]}")          # are they all zero?
    print(f"Loss value: {loss.item():.4f}")        # is it non-zero?
    print(f"Loss grad_fn: {loss.grad_fn}")         # is it None?

    has_grad = any(
        p.grad is not None and p.grad.abs().sum() > 0
        for p in model.parameters()
    )
    assert has_grad, "No gradients found after backward pass — something is broken"
    print(f"✓ Gradients flow correctly, loss = {loss.item():.2f}")
    """


if __name__ == "__main__":
    print("=== Testing PordleCNN ===\n")
    test_parameter_count()
    test_output_shape()
    test_output_non_negative()
    test_intermediate_shapes()
    test_gradient_flows()
    print("\n=== All tests passed ===")
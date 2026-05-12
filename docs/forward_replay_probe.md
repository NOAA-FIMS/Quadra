# Forward replay probe

This diagnostic tests the missing reusable-tape capability:

```text
x.val mutation
-> forward replay
-> updated primal values
-> updated reverse sweep
```

Expected current behavior:

| quantity | expected |
|---|---|
| y value after mutation | stale |
| gradient after mutation | stale or partially stale |

This probe should fail conceptually until had_quadra supports replayable
graph operations or explicit forward recomputation.

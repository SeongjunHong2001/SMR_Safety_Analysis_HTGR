# 가정 입력 파일 (`config/assumptions/*.yaml`)

`config/geometry.yaml`·`operation.yaml`에서 값이 `null`(미확보)인 입력을 계산을 돌리기 위해 임시로 채우는 파일이다.
원문 값이 아니므로 기본 실행에는 쓰지 않는다.

```bash
python main.py --scenario steady_state --condition system_loop --assumptions config/assumptions/<파일>.yaml
```

## 규칙 (`system/builder.py` `load_inputs`)

- 최상위 키는 `meta`, `geometry`, `operation`만 쓴다. `geometry`·`operation` 아래는 각 config와 같은 구조로 깊게 덮어쓴다.
- config에서 `null`이거나 없는 키만 채운다. 원문 값이 있는 키를 다른 값으로 바꾸려 하면 오류다.
- 적용한 키·값·파일은 실행 결과 메타데이터(`*_meta.json`)의 `assumptions`에 남는다.
  계산 영역 밖 장치에 해당하는 항목은 `assumptions_outside_domain`으로 따로 적는다.
- 각 값에는 `[가정]` 주석과 근거를 단다. 원문 자료가 확보되면 config에 넣고 여기서 지운다.
- `meta.kind: calibration`인 파일은 가정이 아니라 보정 결과다 (`config/calibration/`, `--calibration`).

## 이력

- `provisional_s0.yaml` (2026-11-27 작성, 같은 날 삭제)
  - 처음에는 우회 유로 수력직경, SG He 측 수력직경, 노드 110·300·330·340·350 체적, 340 유로를 [가정]으로 채웠다.
  - 수정 기록 #31·#32에서 모두 R19(PSID)·R1 치수로 계산하거나 도면을 해석해 `geometry.yaml`에 넣었다.
  - 남은 항목이 없어 파일을 지웠다. 내용은 커밋 `bd516d1`에 있다.

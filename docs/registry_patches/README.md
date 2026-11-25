# 변수·방정식 변경안 제출 규칙

각 역할은 `roleN/variable_patch.csv`와 `roleN/equation_patch.csv`에
자기 변경안을 기록한다. Role 5만 `docs/variable_registry.csv`,
`docs/equation_registry.csv` 및 Excel master를 통합한다. 기존 이름과
의미를 먼저 검색하고, 중복이면 `update`로 기록한다.

`variable_patch.csv`에는 `action`(add/update), 기존/신규 분류,
SI 단위, owner, Registry 필요 여부, 런타임 키 패턴, 코드 참조,
출처와 사유를 쓴다. 형상·초깃값·실제 수치는 증거가 있을 때만
master에 반영한다. `equation_patch.csv`에는 component/equation_id,
식, 이산화, 미지수, 의존성, 잔차 수, `s_R`, 경계, 근거와 구현 상태를
적는다. 종속된 보존식은 독립 잔차로 중복 등록하지 않는다.

빈 템플릿은 같은 폴더의 `variable_patch_template.csv`,
`equation_patch_template.csv`다. Role 5의 이번 변경 원장은
`role5/`에 있다. 신규 11개 모두 값이 미정이며 실제 원자로의
물리 결과를 뜻하지 않는다.
